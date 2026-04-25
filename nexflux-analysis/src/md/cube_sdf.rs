//! 团簇空间分布函数（Cluster SDF）— mol-cube 第四模式
//!
//! 从 MD 轨迹中识别指定 Qn 等级的磷酸盐团簇（按连通分量的最高个人 Qn 划分），
//! 对每个团簇做刚体对齐（Kabsch + P 原子排列枚举），在局部笛卡尔坐标系叠加统计
//! 各原子类型的三维概率密度，按签名分组输出 Gaussian cube 格式。
//!
//! **签名**：同一帧内团簇的原子类型组合（如 `"Ob:3|On:6|P1:3|P3:1|Zn:2"`），
//! 相同签名归为同一 family，各自取第一个团簇为参考结构。
//!
//! **对齐策略**：
//! - Q0（单 P）：枚举所有 O 原子排列（≤4! = 24），Kabsch 对齐 O 原子
//! - Q1/Q2/Q3：枚举同标签 P 原子组内的所有排列（≤3! = 6），Kabsch 对齐 P 原子
//!
//! **RMSD 质检**：对齐后计算 P 原子（Q0 时为所有 O）的均方根偏差；
//! 超过 `rmsd_warn_threshold` 时打印警告，提示该团簇结构异常。
//!
//! CLI 参数（供 nexflux-cli 接入）：
//!   `--qn 3`             — 目标 Qn 等级 (0/1/2/3)
//!   `--former P`         — 网络形成体元素
//!   `--ligand O`         — 配体元素
//!   `--cutoff-fl 2.4`    — former-ligand 截断半径 [Å]
//!   `--modifier Zn`      — 修饰体元素（可省略）
//!   `--cutoff-ml 2.8`    — modifier-ligand 截断半径 [Å]
//!   `--grid-res 0.1`     — 网格分辨率 [Å/voxel]
//!   `--sigma 1.5`        — 高斯展宽 [voxels]
//!   `--padding 3.0`      — 网格边界留白 [Å]
//!   `--rmsd-warn 0.5`    — RMSD 警告阈值 [Å]

use std::collections::{BTreeMap, HashMap, HashSet};

use nexflux_core::{Atom, Cell, CubeData, Frame, Trajectory};
use nalgebra::{Matrix3, Vector3};
use ndarray::Array3;
use petgraph::{
    graph::{NodeIndex, UnGraph},
    visit::Bfs,
};

// ─── 公开参数 ─────────────────────────────────────────────────────────────────

/// 团簇 SDF 计算参数。
#[derive(Debug, Clone)]
pub struct ClusterSdfParams {
    /// 网络形成体元素，如 `"P"`
    pub former: String,
    /// 配体元素，如 `"O"`
    pub ligand: String,
    /// 目标团簇 Qn 等级（0/1/2/3）；以连通分量内最高个人 Qn 判定
    pub target_qn: u8,
    /// former-ligand 截断半径 [Å]
    pub former_ligand_cutoff: f64,
    /// 修饰体元素（如 `Some("Zn")`），`None` 则不纳入团簇
    pub modifier: Option<String>,
    /// modifier-ligand 截断半径 [Å]（`modifier` 为 `None` 时忽略）
    pub modifier_cutoff: f64,
    /// 网格分辨率 [Å/voxel]
    pub grid_res: f64,
    /// 高斯展宽 sigma [voxels]
    pub sigma: f64,
    /// 网格边界额外留白 [Å]
    pub padding: f64,
    /// RMSD 超过此值时打印警告 [Å]
    pub rmsd_warn_threshold: f64,
}

impl Default for ClusterSdfParams {
    fn default() -> Self {
        Self {
            former: "P".into(),
            ligand: "O".into(),
            target_qn: 3,
            former_ligand_cutoff: 2.4,
            modifier: None,
            modifier_cutoff: 2.8,
            grid_res: 0.1,
            sigma: 1.5,
            padding: 3.0,
            rmsd_warn_threshold: 0.5,
        }
    }
}

// ─── 公开结果 ─────────────────────────────────────────────────────────────────

/// 对齐质量统计（仅含参与对齐的骨架原子的 RMSD）。
#[derive(Debug, Clone, Default)]
pub struct RmsdStats {
    pub mean: f64,
    pub max: f64,
    /// 超出 `rmsd_warn_threshold` 的团簇数
    pub n_warned: usize,
}

/// 同签名（相同原子类型组合）团簇的 SDF 结果。
#[derive(Debug)]
pub struct ClusterFamily {
    /// 签名字符串，如 `"Ob:3|On:6|P1:3|P3:1|Zn:2"`
    pub signature: String,
    /// 各原子类型标签 → 三维密度格点（与所有 family 共享同一网格 origin/spacing）
    pub grids: HashMap<String, CubeData>,
    /// 纳入统计的团簇数量（含参考结构自身）
    pub n_clusters: usize,
    /// 对齐 RMSD 统计（参考结构自身 RMSD=0 计入均值）
    pub rmsd_stats: RmsdStats,
}

/// [`calc_cluster_sdf`] 的返回值。
pub struct ClusterSdfResult {
    /// 签名 → ClusterFamily；可能存在多种签名（拓扑异构体）
    pub families: HashMap<String, ClusterFamily>,
    /// 处理的总帧数（有 Cell 的帧）
    pub n_frames: usize,
    /// 所有 family 的团簇总数
    pub n_clusters_total: usize,
}

// ─── 内部类型 ─────────────────────────────────────────────────────────────────

/// 已提取并在局部坐标系（锚原子为原点）表示的单个团簇。
#[derive(Clone)]
struct ClusterSnapshot {
    /// 各原子类型标签：`"P0"/"P1"/"P2"/"P3"/"Of"/"On"/"Ob"/"Zn"`
    types: Vec<String>,
    /// 局部笛卡尔坐标，锚原子为原点 [Å]
    positions: Vec<Vector3<f64>>,
    /// 锚原子在 `types`/`positions` 中的下标（保留供后续 CLI 使用）
    #[allow(dead_code)]
    anchor_idx: usize,
}

impl ClusterSnapshot {
    /// 按 BTreeMap 排序的原子类型计数字符串，用作 HashMap key。
    fn signature(&self) -> String {
        let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
        for t in &self.types {
            *counts.entry(t.as_str()).or_default() += 1;
        }
        counts.iter()
            .map(|(t, c)| format!("{t}:{c}"))
            .collect::<Vec<_>>()
            .join("|")
    }

    /// 所有 P 原子在 positions 中的下标。
    fn p_indices(&self) -> Vec<usize> {
        self.types.iter().enumerate()
            .filter(|(_, t)| t.starts_with('P'))
            .map(|(i, _)| i)
            .collect()
    }

    /// 所有非修饰体（非 Zn）的 O 原子下标。
    fn o_indices(&self) -> Vec<usize> {
        self.types.iter().enumerate()
            .filter(|(_, t)| matches!(t.as_str(), "Of" | "On" | "Ob"))
            .map(|(i, _)| i)
            .collect()
    }
}

/// 按签名累积的统计数据。
struct FamilyAccumulator {
    reference: ClusterSnapshot,
    /// type_label → 所有对齐后的原子坐标（跨所有帧和团簇）
    positions_by_type: HashMap<String, Vec<Vector3<f64>>>,
    n_clusters: usize,
    rmsd_sum: f64,
    rmsd_max: f64,
    n_warned: usize,
}

impl FamilyAccumulator {
    /// 以第一个团簇为参考结构初始化；参考自身以 RMSD=0 计入统计。
    fn new(reference: ClusterSnapshot) -> Self {
        let mut acc = Self {
            positions_by_type: HashMap::new(),
            reference: reference.clone(),
            n_clusters: 0,
            rmsd_sum: 0.0,
            rmsd_max: 0.0,
            n_warned: 0,
        };
        acc.push(&reference, 0.0, f64::MAX); // 参考本身不触发警告
        acc
    }

    fn push(&mut self, snapshot: &ClusterSnapshot, rmsd: f64, threshold: f64) {
        for (t, &p) in snapshot.types.iter().zip(&snapshot.positions) {
            self.positions_by_type.entry(t.clone()).or_default().push(p);
        }
        self.n_clusters += 1;
        self.rmsd_sum += rmsd;
        if rmsd > self.rmsd_max { self.rmsd_max = rmsd; }
        if rmsd > threshold {
            self.n_warned += 1;
            eprintln!(
                "[ClusterSDF] 警告：RMSD = {rmsd:.3} Å > 阈值 {threshold:.3} Å，\
                 团簇结构偏差过大，请人工甄别。"
            );
        }
    }
}

// ─── 主入口 ───────────────────────────────────────────────────────────────────

/// 计算指定 Qn 团簇的空间分布函数（SDF）。
///
/// 返回 `None` 若：
/// - 轨迹无帧或每帧无 Cell
/// - 整条轨迹未找到任何目标 Qn 团簇
pub fn calc_cluster_sdf(
    traj: &Trajectory,
    params: &ClusterSdfParams,
) -> Option<ClusterSdfResult> {
    if traj.frames.is_empty() { return None; }

    let mut accumulators: HashMap<String, FamilyAccumulator> = HashMap::new();
    let mut n_frames = 0usize;

    for frame in &traj.frames {
        let cell = match frame.cell.as_ref() {
            Some(c) => c,
            None => continue,
        };
        n_frames += 1;

        for mut snapshot in process_frame(frame, cell, params) {
            let sig = snapshot.signature();
            if let Some(acc) = accumulators.get_mut(&sig) {
                let rmsd = align_to_reference(&mut snapshot, &acc.reference, params.target_qn);
                acc.push(&snapshot, rmsd, params.rmsd_warn_threshold);
            } else {
                accumulators.insert(sig, FamilyAccumulator::new(snapshot));
            }
        }
    }

    if accumulators.is_empty() { return None; }

    let n_clusters_total = accumulators.values().map(|a| a.n_clusters).sum();
    let families = accumulators
        .into_iter()
        .map(|(sig, acc)| {
            let family = build_family(acc, params);
            (sig, family)
        })
        .collect();

    Some(ClusterSdfResult { families, n_frames, n_clusters_total })
}

// ─── 逐帧处理 ─────────────────────────────────────────────────────────────────

fn process_frame(frame: &Frame, cell: &Cell, params: &ClusterSdfParams) -> Vec<ClusterSnapshot> {
    // 1. 按元素收集全局原子索引
    let mut p_global: Vec<usize> = Vec::new();
    let mut o_global: Vec<usize> = Vec::new();
    let mut mod_global: Vec<usize> = Vec::new();

    for (i, atom) in frame.atoms.iter().enumerate() {
        if atom.element == params.former          { p_global.push(i); }
        else if atom.element == params.ligand     { o_global.push(i); }
        else if params.modifier.as_deref() == Some(atom.element.as_str()) {
            mod_global.push(i);
        }
    }

    let np = p_global.len();
    let no = o_global.len();
    if np == 0 || no == 0 { return Vec::new(); }

    // 2. P-O 邻接表（局部下标）
    let fl_cutoff2 = params.former_ligand_cutoff * params.former_ligand_cutoff;
    let mut p_o_adj: Vec<Vec<usize>> = vec![Vec::new(); np]; // p_local → [o_local]
    let mut o_p_adj: Vec<Vec<usize>> = vec![Vec::new(); no]; // o_local → [p_local]

    for (pi, &pa_idx) in p_global.iter().enumerate() {
        let pa_pos = frame.atoms[pa_idx].position;
        for (oi, &oa_idx) in o_global.iter().enumerate() {
            let diff = cell.minimum_image(frame.atoms[oa_idx].position - pa_pos);
            if diff.norm_squared() < fl_cutoff2 {
                p_o_adj[pi].push(oi);
                o_p_adj[oi].push(pi);
            }
        }
    }

    // 3. O 原子类型：Of=0P, On=1P, Ob≥2P
    let o_type: Vec<&str> = o_p_adj.iter().map(|ps| match ps.len() {
        0 => "Of",
        1 => "On",
        _ => "Ob",
    }).collect();

    // 4. 每个 P 的个人 Qn（连接的 Ob 数量）
    let p_qn: Vec<u8> = (0..np)
        .map(|pi| p_o_adj[pi].iter().filter(|&&oi| o_type[oi] == "Ob").count() as u8)
        .collect();

    // 5. 建 P-P 图（共享 Ob 则有边），BFS 求连通分量
    let mut g: UnGraph<(), ()> = UnGraph::with_capacity(np, np * 2);
    let nodes: Vec<NodeIndex> = (0..np).map(|_| g.add_node(())).collect();

    for oi in 0..no {
        if o_type[oi] == "Ob" {
            let ps = &o_p_adj[oi];
            for i in 0..ps.len() {
                for j in (i + 1)..ps.len() {
                    g.add_edge(nodes[ps[i]], nodes[ps[j]], ());
                }
            }
        }
    }

    let mut visited = vec![false; np];
    let mut components: Vec<Vec<usize>> = Vec::new(); // p_local 下标组

    for start in 0..np {
        if visited[start] { continue; }
        let mut comp: Vec<usize> = Vec::new();
        let mut bfs = Bfs::new(&g, nodes[start]);
        while let Some(nx) = bfs.next(&g) {
            let pi = nx.index();
            if !visited[pi] {
                visited[pi] = true;
                comp.push(pi);
            }
        }
        components.push(comp);
    }

    // 6. 筛选 target_qn 团簇，提取快照
    let ml_cutoff2 = params.modifier_cutoff * params.modifier_cutoff;

    components.into_iter().filter_map(|comp_p| {
        let max_qn = comp_p.iter().map(|&pi| p_qn[pi]).max()?;
        if max_qn != params.target_qn { return None; }

        extract_snapshot(
            frame, cell,
            &comp_p, &p_global, &o_global, &mod_global,
            &p_o_adj, &o_type, &p_qn,
            ml_cutoff2, params.target_qn,
        )
    }).collect()
}

/// 从单个连通分量提取团簇快照（局部坐标，锚原子为原点）。
fn extract_snapshot(
    frame: &Frame,
    cell: &Cell,
    comp_p: &[usize],
    p_global: &[usize],
    o_global: &[usize],
    mod_global: &[usize],
    p_o_adj: &[Vec<usize>],
    o_type: &[&str],
    p_qn: &[u8],
    ml_cutoff2: f64,
    target_qn: u8,
) -> Option<ClusterSnapshot> {
    // 收集该团簇的 O 原子（连接到任一 P 的 O）
    let mut cluster_o: Vec<usize> = {
        let mut set: HashSet<usize> = HashSet::new();
        for &pi in comp_p {
            set.extend(p_o_adj[pi].iter().copied());
        }
        let mut v: Vec<usize> = set.into_iter().collect();
        v.sort_unstable();
        v
    };

    // 确定锚原子全局索引与类型标签
    let anchor_ga: usize = if target_qn == 1 {
        // Q1：以唯一 Ob 为锚
        let ob_oi = cluster_o.iter().find(|&&oi| o_type[oi] == "Ob")?;
        o_global[*ob_oi]
    } else {
        // Q0/Q2/Q3：以最高个人 Qn 的 P 为锚
        let anchor_pi = *comp_p.iter().max_by_key(|&&pi| p_qn[pi])?;
        p_global[anchor_pi]
    };

    let anchor_pos = frame.atoms[anchor_ga].position;

    // 收集修饰体原子（紧邻团簇 O 原子）
    let cluster_mod: Vec<usize> = if ml_cutoff2 > 0.0 && !mod_global.is_empty() {
        let cluster_o_pos: Vec<Vector3<f64>> = cluster_o.iter()
            .map(|&oi| frame.atoms[o_global[oi]].position)
            .collect();

        mod_global.iter().enumerate().filter_map(|(mi, &ma_idx)| {
            let ma_pos = frame.atoms[ma_idx].position;
            let is_near = cluster_o_pos.iter().any(|&op| {
                cell.minimum_image(op - ma_pos).norm_squared() < ml_cutoff2
            });
            if is_near { Some(mi) } else { None }
        }).collect()
    } else {
        Vec::new()
    };

    // 构建原子列表：P（锚优先，余按 Qn 降序）→ O（Ob > On > Of）→ modifier
    let mut types: Vec<String> = Vec::new();
    let mut raw_pos: Vec<Vector3<f64>> = Vec::new();
    let mut anchor_idx = 0usize;

    // P 原子
    let mut sorted_p = comp_p.to_vec();
    sorted_p.sort_by(|&a, &b| p_qn[b].cmp(&p_qn[a]).then(a.cmp(&b)));
    for &pi in &sorted_p {
        let ga = p_global[pi];
        if ga == anchor_ga { anchor_idx = types.len(); }
        types.push(format!("P{}", p_qn[pi]));
        raw_pos.push(frame.atoms[ga].position);
    }

    // O 原子
    cluster_o.sort_by(|&a, &b| {
        let ord = |t: &str| match t { "Ob" => 0u8, "On" => 1, _ => 2 };
        ord(o_type[a]).cmp(&ord(o_type[b])).then(a.cmp(&b))
    });
    for &oi in &cluster_o {
        let ga = o_global[oi];
        if ga == anchor_ga { anchor_idx = types.len(); }
        types.push(o_type[oi].to_string());
        raw_pos.push(frame.atoms[ga].position);
    }

    // modifier 原子
    for mi in cluster_mod {
        let ga = mod_global[mi];
        types.push(frame.atoms[ga].element.clone());
        raw_pos.push(frame.atoms[ga].position);
    }

    if types.is_empty() { return None; }

    // PBC 展开：以锚原子为原点，用最小镜像得到局部坐标
    let positions: Vec<Vector3<f64>> = raw_pos.iter()
        .map(|&pos| cell.minimum_image(pos - anchor_pos))
        .collect();

    Some(ClusterSnapshot { types, positions, anchor_idx })
}

// ─── 刚体对齐 ─────────────────────────────────────────────────────────────────

/// 将 `mobile` 对齐到 `reference`，原地旋转 `mobile.positions`，返回骨架 RMSD [Å]。
///
/// - Q0：以 O 原子枚举排列找最优旋转（P 已在原点，自动对齐）
/// - Q1/Q2/Q3：以 P 原子枚举同标签组内排列找最优旋转
fn align_to_reference(
    mobile: &mut ClusterSnapshot,
    reference: &ClusterSnapshot,
    target_qn: u8,
) -> f64 {
    let rot = if target_qn == 0 {
        best_rotation_by_permutation(reference, mobile, mobile.o_indices(), reference.o_indices())
    } else {
        best_rotation_by_permutation(reference, mobile, mobile.p_indices(), reference.p_indices())
    };

    for pos in mobile.positions.iter_mut() {
        *pos = rot * *pos;
    }

    // 骨架 RMSD：Q0 用 O，其余用 P
    let (ref_idx, mob_idx) = if target_qn == 0 {
        (reference.o_indices(), mobile.o_indices())
    } else {
        (reference.p_indices(), mobile.p_indices())
    };

    if ref_idx.is_empty() || ref_idx.len() != mob_idx.len() {
        return 0.0;
    }

    let ref_pts: Vec<_> = ref_idx.iter().map(|&i| reference.positions[i]).collect();
    let mob_pts: Vec<_> = mob_idx.iter().map(|&i| mobile.positions[i]).collect();
    compute_rmsd(&ref_pts, &mob_pts)
}

/// 枚举 `mob_atom_indices` 对应原子的所有排列，对每种排列计算 Kabsch 旋转，
/// 返回使 RMSD 最小的旋转矩阵。
fn best_rotation_by_permutation(
    reference: &ClusterSnapshot,
    mobile: &ClusterSnapshot,
    mob_atom_indices: Vec<usize>,
    ref_atom_indices: Vec<usize>,
) -> Matrix3<f64> {
    if mob_atom_indices.is_empty() || mob_atom_indices.len() != ref_atom_indices.len() {
        return Matrix3::identity();
    }

    let ref_pts: Vec<Vector3<f64>> = ref_atom_indices.iter()
        .map(|&i| reference.positions[i])
        .collect();
    let mob_pool: Vec<Vector3<f64>> = mob_atom_indices.iter()
        .map(|&i| mobile.positions[i])
        .collect();

    let n = mob_pool.len();
    let mut perm: Vec<usize> = (0..n).collect();
    let mut best_rmsd = f64::MAX;
    let mut best_rot = Matrix3::identity();

    heap_permutations(&mut perm, n, &mut |p: &[usize]| {
        let mob_pts: Vec<Vector3<f64>> = p.iter().map(|&i| mob_pool[i]).collect();
        let rot = kabsch_rotation(&ref_pts, &mob_pts);
        let rotated: Vec<_> = mob_pts.iter().map(|&m| rot * m).collect();
        let rmsd = compute_rmsd(&ref_pts, &rotated);
        if rmsd < best_rmsd {
            best_rmsd = rmsd;
            best_rot = rot;
        }
    });

    best_rot
}

/// Kabsch 算法：返回旋转矩阵 R，使 R * mob_i ≈ ref_i（列向量约定）。
///
/// H = Σ mob_i · ref_iᵀ → SVD(H) = U S Vᵀ → R = U · diag(1,1,d) · Vᵀ
/// 其中 d = sign(det(V · Uᵀ))，确保 det(R)=+1（纯旋转，无反射）。
fn kabsch_rotation(ref_pts: &[Vector3<f64>], mob_pts: &[Vector3<f64>]) -> Matrix3<f64> {
    let mut h = Matrix3::zeros();
    for (r, m) in ref_pts.iter().zip(mob_pts) {
        h += m * r.transpose();
    }
    let svd = nalgebra::linalg::SVD::new(h, true, true);
    let (Some(u), Some(vt)) = (svd.u, svd.v_t) else { return Matrix3::identity() };
    let v = vt.transpose();
    let d = (v * u.transpose()).determinant().signum();
    let diag = Matrix3::from_diagonal(&Vector3::new(1.0, 1.0, if d >= 0.0 { 1.0 } else { -1.0 }));
    u * diag * vt
}

/// Heap's algorithm（非递归风格）：枚举 `arr[..k]` 的所有 k! 排列，每次调用 `cb`。
fn heap_permutations(arr: &mut [usize], k: usize, cb: &mut impl FnMut(&[usize])) {
    if k == 1 {
        cb(arr);
        return;
    }
    for i in 0..k {
        heap_permutations(arr, k - 1, cb);
        if k % 2 == 0 {
            arr.swap(i, k - 1);
        } else {
            arr.swap(0, k - 1);
        }
    }
}

fn compute_rmsd(a: &[Vector3<f64>], b: &[Vector3<f64>]) -> f64 {
    if a.is_empty() { return 0.0; }
    let msd: f64 = a.iter().zip(b).map(|(x, y)| (x - y).norm_squared()).sum::<f64>() / a.len() as f64;
    msd.sqrt()
}

// ─── 网格构建 ─────────────────────────────────────────────────────────────────

/// 将 `FamilyAccumulator` 转换为 `ClusterFamily`（构建 3D 直方图 + 高斯展宽）。
fn build_family(acc: FamilyAccumulator, params: &ClusterSdfParams) -> ClusterFamily {
    let sig = acc.reference.signature();

    // 确定全局网格范围（所有类型的所有坐标共用同一 box）
    let all_pos: Vec<Vector3<f64>> = acc.positions_by_type.values()
        .flat_map(|v| v.iter().copied())
        .collect();

    if all_pos.is_empty() {
        return ClusterFamily {
            signature: sig,
            grids: HashMap::new(),
            n_clusters: acc.n_clusters,
            rmsd_stats: finalize_rmsd(&acc),
        };
    }

    let mut min_xyz = all_pos[0];
    let mut max_xyz = all_pos[0];
    for p in &all_pos {
        min_xyz.x = min_xyz.x.min(p.x);
        min_xyz.y = min_xyz.y.min(p.y);
        min_xyz.z = min_xyz.z.min(p.z);
        max_xyz.x = max_xyz.x.max(p.x);
        max_xyz.y = max_xyz.y.max(p.y);
        max_xyz.z = max_xyz.z.max(p.z);
    }

    let origin = min_xyz - Vector3::repeat(params.padding);
    let box_size = (max_xyz - min_xyz) + Vector3::repeat(2.0 * params.padding);
    let nx = (box_size.x / params.grid_res).ceil() as usize + 1;
    let ny = (box_size.y / params.grid_res).ceil() as usize + 1;
    let nz = (box_size.z / params.grid_res).ceil() as usize + 1;

    let spacing = Matrix3::from_diagonal(&Vector3::repeat(params.grid_res));
    let ref_frame = build_reference_frame(&acc.reference, origin, box_size);

    // 每种 type_label 独立建格点，共享 origin/spacing/ref_frame
    let rmsd_stats = finalize_rmsd(&acc);
    let n_clusters = acc.n_clusters;

    let grids: HashMap<String, CubeData> = acc.positions_by_type.iter().map(|(label, positions)| {
        let cube = build_cube_for_type(
            positions, &origin, nx, ny, nz, params.grid_res, params.sigma,
            spacing, ref_frame.clone(),
        );
        (label.clone(), cube)
    }).collect();

    ClusterFamily { signature: sig, grids, n_clusters, rmsd_stats }
}

/// 对单个原子类型建 3D 密度格点（直方图 + 高斯展宽 + 均值归一化）。
fn build_cube_for_type(
    positions: &[Vector3<f64>],
    origin: &Vector3<f64>,
    nx: usize,
    ny: usize,
    nz: usize,
    grid_res: f64,
    sigma: f64,
    spacing: Matrix3<f64>,
    frame: Frame,
) -> CubeData {
    let mut counts = Array3::<f64>::zeros((nx, ny, nz));

    for &pos in positions {
        let rel = pos - origin;
        let ix = (rel.x / grid_res).floor() as isize;
        let iy = (rel.y / grid_res).floor() as isize;
        let iz = (rel.z / grid_res).floor() as isize;
        if ix >= 0 && iy >= 0 && iz >= 0
            && (ix as usize) < nx && (iy as usize) < ny && (iz as usize) < nz
        {
            counts[[ix as usize, iy as usize, iz as usize]] += 1.0;
        }
    }

    let smoothed = gaussian_filter_3d(counts, sigma);

    // 均值归一化（mean of all voxels = 1.0），便于可视化对比
    let mean = smoothed.mean().unwrap_or(1.0);
    let data = if mean > 0.0 { smoothed.mapv(|v| v / mean) } else { smoothed };

    CubeData { frame, data, origin: *origin, spacing }
}

/// 构建放置在 cube box 中心的参考结构帧（用于 cube 文件头）。
fn build_reference_frame(
    snapshot: &ClusterSnapshot,
    origin: Vector3<f64>,
    box_size: Vector3<f64>,
) -> Frame {
    let center = origin + box_size * 0.5;
    let mut frame = Frame::new();
    for (t, &pos) in snapshot.types.iter().zip(&snapshot.positions) {
        let elem = type_to_element(t);
        let shifted = pos + center; // anchor 在原点，平移到 box 中心
        frame.add_atom(Atom::new(elem, shifted));
    }
    frame
}

/// 将内部 type label 转为化学元素符号。
fn type_to_element(label: &str) -> &str {
    if label.starts_with('P') { "P" }
    else if matches!(label, "Of" | "On" | "Ob") { "O" }
    else { label } // modifier 直接用元素符号（如 "Zn"）
}

fn finalize_rmsd(acc: &FamilyAccumulator) -> RmsdStats {
    let mean = if acc.n_clusters > 0 { acc.rmsd_sum / acc.n_clusters as f64 } else { 0.0 };
    RmsdStats { mean, max: acc.rmsd_max, n_warned: acc.n_warned }
}

// ─── 高斯展宽 ─────────────────────────────────────────────────────────────────

/// 三维高斯滤波（可分离，沿三轴依次做 1D 卷积）。
fn gaussian_filter_3d(data: Array3<f64>, sigma: f64) -> Array3<f64> {
    if sigma <= 0.0 { return data; }
    let kernel = gaussian_kernel(sigma);
    let out = convolve1d_axis(&data, &kernel, 0);
    let out = convolve1d_axis(&out, &kernel, 1);
    convolve1d_axis(&out, &kernel, 2)
}

/// 归一化 1D 高斯核（截断于 ±4σ）。
fn gaussian_kernel(sigma: f64) -> Vec<f64> {
    let radius = (4.0 * sigma).ceil() as isize;
    let mut kernel: Vec<f64> = (-radius..=radius)
        .map(|i| (-(i * i) as f64 / (2.0 * sigma * sigma)).exp())
        .collect();
    let sum: f64 = kernel.iter().sum();
    kernel.iter_mut().for_each(|v| *v /= sum);
    kernel
}

/// 沿指定轴（0/1/2）做 zero-padding 1D 卷积。
fn convolve1d_axis(data: &Array3<f64>, kernel: &[f64], axis: usize) -> Array3<f64> {
    let (nx, ny, nz) = (data.dim().0, data.dim().1, data.dim().2);
    let mut out = Array3::<f64>::zeros((nx, ny, nz));
    let half = (kernel.len() / 2) as isize;
    let len = [nx, ny, nz][axis] as isize;

    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                let center = [ix as isize, iy as isize, iz as isize][axis];
                let val: f64 = kernel.iter().enumerate().map(|(ki, &kv)| {
                    let pos = center + (ki as isize) - half;
                    if pos < 0 || pos >= len { return 0.0; }
                    let idx = match axis {
                        0 => [pos as usize, iy, iz],
                        1 => [ix, pos as usize, iz],
                        _ => [ix, iy, pos as usize],
                    };
                    data[idx] * kv
                }).sum();
                out[[ix, iy, iz]] = val;
            }
        }
    }
    out
}

// ─── 测试 ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    #[test]
    fn test_signature_sorted() {
        let snap = ClusterSnapshot {
            types: vec!["P1".into(), "Ob".into(), "P3".into(), "On".into()],
            positions: vec![Vector3::zeros(); 4],
            anchor_idx: 2,
        };
        assert_eq!(snap.signature(), "Ob:1|On:1|P1:1|P3:1");
    }

    #[test]
    fn test_kabsch_identity() {
        let pts: Vec<Vector3<f64>> = vec![
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let rot = kabsch_rotation(&pts, &pts);
        for p in &pts {
            let rp = rot * p;
            assert!((rp - p).norm() < 1e-9, "identity rotation failed");
        }
    }

    #[test]
    fn test_kabsch_known_rotation() {
        // 90° around z-axis: (1,0,0) → (0,1,0)
        let ref_pts = vec![Vector3::new(0.0, 1.0, 0.0), Vector3::new(-1.0, 0.0, 0.0)];
        let mob_pts = vec![Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, -1.0, 0.0)];
        let rot = kabsch_rotation(&ref_pts, &mob_pts);
        for (r, m) in ref_pts.iter().zip(&mob_pts) {
            let aligned = rot * m;
            assert!((aligned - r).norm() < 1e-8, "rotation mismatch: {aligned} vs {r}");
        }
    }

    #[test]
    fn test_kabsch_no_reflection() {
        // 确保 det(R) = +1（纯旋转，非反射）
        let ref_pts = vec![
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let mob_pts = vec![
            Vector3::new(-1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let rot = kabsch_rotation(&ref_pts, &mob_pts);
        assert!((rot.determinant() - 1.0).abs() < 1e-9, "det(R) != 1");
    }

    #[test]
    fn test_compute_rmsd_zero() {
        let pts = vec![Vector3::new(1.0, 2.0, 3.0), Vector3::new(4.0, 5.0, 6.0)];
        assert!((compute_rmsd(&pts, &pts)).abs() < 1e-12);
    }

    #[test]
    fn test_heap_permutations_count() {
        let mut arr: Vec<usize> = (0..3).collect();
        let mut count = 0usize;
        heap_permutations(&mut arr, 3, &mut |_| count += 1);
        assert_eq!(count, 6); // 3! = 6
    }

    #[test]
    fn test_heap_permutations_all_unique() {
        let mut arr: Vec<usize> = (0..3).collect();
        let mut perms: Vec<Vec<usize>> = Vec::new();
        heap_permutations(&mut arr, 3, &mut |p| perms.push(p.to_vec()));
        assert_eq!(perms.len(), 6);
        let unique: HashSet<Vec<usize>> = perms.into_iter().collect();
        assert_eq!(unique.len(), 6);
    }

    #[test]
    fn test_gaussian_kernel_normalized() {
        let k = gaussian_kernel(2.0);
        let sum: f64 = k.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10, "kernel sum = {sum}");
    }

    #[test]
    fn test_gaussian_filter_preserves_nonzero() {
        let mut data = Array3::<f64>::zeros((10, 10, 10));
        data[[5, 5, 5]] = 1.0;
        let out = gaussian_filter_3d(data, 1.0);
        assert!(out[[5, 5, 5]] > 0.0);
        assert!(out[[4, 5, 5]] > 0.0);
        assert!(out[[5, 5, 4]] > 0.0);
    }

    #[test]
    fn test_type_to_element() {
        assert_eq!(type_to_element("P3"), "P");
        assert_eq!(type_to_element("Ob"), "O");
        assert_eq!(type_to_element("On"), "O");
        assert_eq!(type_to_element("Zn"), "Zn");
    }
}
