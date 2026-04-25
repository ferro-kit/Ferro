# molflow 项目架构说明

## 依赖层次

```
molflow-cli / molflow-python        ← 唯一允许组合多 crate 的层
    ├── molflow-core                ← 纯数据结构，无 IO 逻辑
    ├── molflow-io        → core    ← 格式读写
    ├── molflow-analysis  → core    ← 后处理分析
    └── molflow-workflow  → core    ← QC 软件输入文件生成
```

**规则**：中间层 crate（io / analysis / workflow）不得互相依赖。

---

## molflow-core

纯数据结构与静态参考数据，无 IO 逻辑，无外部 IO 依赖。

```
molflow-core/src/
├── lib.rs
├── atom.rs          # Atom { element, position, label, mass, magmom, charge }
├── frame.rs         # Frame { atoms, cell, pbc, energy, forces, velocities, ... }
├── trajectory.rs    # Trajectory { frames, metadata }; tail(), first() 等
├── cell.rs          # Cell(Matrix3<f64>); lengths/angles/volume/wrap/minimum_image
├── units.rs         # LengthUnit/EnergyUnit/PressureUnit/TimeUnit + 转换函数
├── error.rs         # ChemError (thiserror)
└── data/
    ├── elements.rs  # 静态元素表：原子序数、质量、氧化态、电负性
    └── compounds.rs # 常见化合物参考数据（分子量、密度）供 molflow-structure 使用
```

关键设计决策：
- **`Molecule` 类型不存在**，`Frame` 通过 `pbc = [false; 3]` 表示非周期体系
- 原子 index 由位置（`Vec<Atom>` 下标）隐式决定，不存储 `index` 字段
- NPT 轨迹的逐帧盒子变化由 `Frame.cell: Option<Cell>` 自然处理

---

## molflow-io

```
molflow-io/src/
├── lib.rs
├── readers/
│   ├── mod.rs
│   ├── xyz.rs           # read_xyz
│   ├── extxyz.rs        # read_extxyz（支持 Lattice/Properties 头部）
│   ├── pdb.rs           # read_pdb（MODEL/ENDMDL 多帧）
│   ├── cif.rs           # read_cif
│   ├── poscar.rs        # read_poscar / read_contcar
│   ├── lammps_data.rs   # read_lammps_data
│   ├── lammps_dump.rs   # read_lammps_dump（流式逐帧解析）
│   ├── cp2k_inp.rs      # read_cp2k_inp
│   ├── cp2k_restart.rs  # read_cp2k_restart
│   └── qe_input.rs      # read_qe_input
└── writers/
    ├── mod.rs
    ├── xyz.rs / extxyz.rs / pdb.rs / cif.rs / poscar.rs
    ├── lammps_data.rs / lammps_dump.rs
    ├── qe_input.rs
    └── cube.rs          # write_cube（Gaussian cube 格式）
```

所有 reader 返回 `Result<Trajectory>`；所有 writer 接受 `&Trajectory`。

格式自动检测在 `molflow-cli/src/io_dispatch.rs`，按文件名（POSCAR/CONTCAR）和扩展名路由。

---

## molflow-analysis

```
molflow-analysis/src/
├── lib.rs
├── geometry.rs              # 键长、键角、二面角、回转半径、包围盒
├── trajectory_analysis.rs   # RMSD、质心轨迹
├── properties.rs            # 偶极矩等
├── md/
│   ├── mod.rs
│   ├── gr.rs                # g(r) + CN(r) + 键长统计；rayon 并行
│   ├── sq.rs                # S(q)（Faber-Ziman；XRD/neutron 加权）
│   ├── msd.rs               # MSD(t)，各轴分量
│   ├── angle.rs             # 键角分布 P(θ)
│   ├── vanhove.rs           # Van Hove 自相关函数 Gs(r,τ)
│   ├── vacf.rs              # 速度自相关函数 C_v(t) + Green-Kubo 积分
│   ├── rotcorr.rs           # 转动相关函数 C₂(t)（P₂ Legendre）
│   ├── cube_density.rs      # 空间数密度 / 速度 / 力场（cube 格式）
│   └── scattering_data.rs   # X 射线 / 中子散射形状因子表
└── network/
    ├── mod.rs               # NetworkParams, NetworkResult, calc_network
    ├── cn.rs                # build_neighbor_map, collect_cn, collect_cn_total
    ├── ligand_class.rs      # build_ligand_nf_map, classify_ligands (FO/NBO/BO/OBO)
    └── qn.rs                # calc_qn → Qn(mX) 标签
```

### network 模块设计要点

**输入参数**：`NetworkParams { cutoffs: HashMap<(former, ligand), f64> }`

**三层分析**：
```
截断半径邻居搜索 (CN)
        ↓
配体原子分类 (FO / NBO / BO / OBO)
        ↓
Qn 物种分布 (Q^n(mX))
```

**配体分类规则**：
| NF 邻居数 | 标签示例 |
|---|---|
| 0 | `FO` |
| 1 | `NBO(P)` |
| 2 | `BO(P-Zn)`（字母序） |
| ≥3 | `OBO(P-Zn-Zn)`（字母序）|

**Qn(mX) 定义**：
- n = 该 former 原子的所有配体中，属于 BO 或 OBO 类型的数量
- mX = n 个桥接配体中，通过异元素 X（X ≠ former）桥接的数量
- 每个桥接配体对每种异元素只计一次（确保 mX ≤ n）

---

## molflow-workflow

```
molflow-workflow/src/
└── job_builder.rs   # GaussianJobBuilder; 接受 Trajectory，输出 .gjf
```

---

## molflow-cli

```
molflow-cli/src/
├── lib.rs
├── main.rs              # molflow 总览帮助
├── io_dispatch.rs       # read_trajectory / write_trajectory（格式自动检测）
├── help.rs              # mol-traj / mol-corr / mol-cube 模式专属帮助
├── args/
│   ├── mod.rs
│   ├── common.rs        # CommonArgs（-i/-o/--last-n/--ncore）
│   ├── traj.rs          # TrajMode + SqWeightingCli
│   ├── corr.rs          # CorrMode
│   └── cube.rs          # CubeCliMode
└── bin/
    ├── convert.rs       # mol-convert
    ├── info.rs          # mol-info
    ├── job.rs           # mol-job
    ├── traj.rs          # mol-traj -m gr|sq|msd|angle
    ├── corr.rs          # mol-corr -m vacf|rotcorr|vanhove
    ├── cube.rs          # mol-cube -m density|velocity|force
    └── network.rs       # mol-network --Former-Ligand=cutoff
```

### mol-network 参数解析

`--P-O=2.5` 等 pair 参数首字母大写，与 clap 普通参数（均为小写）在预处理阶段分离，
再分别交给自定义解析器和 clap 处理，避免 clap 拒绝未知参数名。

---

## molflow-python（暂未启用）

```
molflow-python/src/
├── lib.rs      # #[pymodule] 入口
├── types.rs    # PyTrajectory（#[pyclass] 包裹 Trajectory）
├── io.rs / structure.rs / analysis.rs
```

零 Python 感知的纯 Rust 库；所有 PyO3 胶水代码只在此 crate。

---

## 依赖清单

| Crate | 用途 |
|---|---|
| `nalgebra` | 坐标 Vector3<f64>，晶胞 Matrix3<f64> |
| `ndarray` | 多维数组（体素网格等） |
| `rayon` | 数据并行（gr、cube 等） |
| `thiserror` | ChemError derive |
| `anyhow` | CLI 错误传播 |
| `clap` | CLI 参数解析 |
| `rust_xlsxwriter` | mol-network xlsx 输出 |
| `pyo3` | Python 绑定（molflow-python 专用） |

---

## 扩展指南

### 添加文件格式
1. `molflow-io/src/readers/<fmt>.rs` — 返回 `Result<Trajectory>`
2. `molflow-io/src/writers/<fmt>.rs`
3. 在 `readers/mod.rs`、`writers/mod.rs` 中导出
4. 在 `molflow-cli/src/io_dispatch.rs` 中添加扩展名路由

### 添加分析功能
1. 在 `molflow-analysis/src/` 实现（接受/返回 `Trajectory` 或中间结果）
2. 在 `molflow-analysis/src/lib.rs` re-export
3. 在对应 CLI bin 中接入

### 添加 QC 软件支持
1. `molflow-workflow/src/job_builder.rs` 实现 builder
2. 在 `molflow-cli/src/bin/job.rs` 添加软件分支
