# molflow

用 Rust 编写的模块化计算化学工具包，主要面向**周期性体系**（晶体、表面、玻璃）。

## 项目结构

```
molflow/
├── molflow-core/       # 核心数据结构（Atom, Frame, Trajectory, Cell）
├── molflow-io/         # 文件读写（10 种格式）
├── molflow-analysis/   # 后处理分析（MD 分析、玻璃网络结构分析）
├── molflow-workflow/   # 计算输入文件生成（Gaussian 等）
├── molflow-cli/        # 命令行工具集
└── molflow-python/     # Python 绑定（PyO3，暂未启用）
```

依赖层次（严格单向，中间层不得互相依赖）：
```
molflow-cli / molflow-python
    ├── molflow-core
    ├── molflow-io        → core
    ├── molflow-analysis  → core
    └── molflow-workflow  → core
```

## 安装与编译

```bash
cargo build --release
```

编译产物在 `target/release/`，包含以下可执行文件：
`mol-convert`、`mol-info`、`mol-job`、`mol-traj`、`mol-corr`、`mol-cube`、`mol-network`

## 支持的文件格式

| 格式 | 读 | 写 |
|---|:---:|:---:|
| XYZ | ✓ | ✓ |
| Extended XYZ | ✓ | ✓ |
| PDB | ✓ | ✓ |
| CIF | ✓ | ✓ |
| POSCAR / CONTCAR | ✓ | ✓ |
| LAMMPS data | ✓ | ✓ |
| LAMMPS dump | ✓ | ✓ |
| CP2K .inp | ✓ | — |
| CP2K .restart | ✓ | — |
| Quantum ESPRESSO .in | ✓ | ✓ |
| Gaussian cube | — | ✓ |

## 命令行工具

### mol-convert — 格式转换
```bash
mol-convert -i structure.xyz -o structure.pdb
mol-convert -i traj.dump -o traj.extxyz
```

### mol-info — 结构信息
```bash
mol-info -i water.xyz
mol-info -i NaCl.cif
```

### mol-job — 生成计算输入文件
```bash
mol-job -i water.xyz -s gaussian -m B3LYP -b "6-31G*" -o job.gjf
```

### mol-traj — 轨迹结构分析
```bash
mol-traj -m gr    -i traj.dump                          # 径向分布函数 g(r)
mol-traj -m sq    -i traj.dump --weighting both         # 结构因子 S(q)
mol-traj -m msd   -i traj.dump --dt 2.0 --elements Li   # 均方位移
mol-traj -m angle -i traj.dump --r-cut-ab 2.0           # 键角分布
mol-traj -m gr                                          # 无 -i → 显示详细帮助
```

### mol-corr — 相关函数
```bash
mol-corr -m vacf    -i traj.dump --dt 2.0              # 速度自相关函数
mol-corr -m rotcorr -i traj.xyz --center O --neighbor H # 转动相关函数
mol-corr -m vanhove -i traj.dump --tau 100              # Van Hove 自相关函数
```

### mol-cube — 空间密度分布
```bash
mol-cube -m density  -i traj.dump                      # 数密度 (原子/Å³)
mol-cube -m velocity -i traj.dump                      # 速度场
mol-cube -m force    -i traj.dump --elements O         # 力场
```
输出为 Gaussian cube 格式，可直接用 VESTA / VMD 可视化。

### mol-network — 玻璃网络结构分析
针对周期性体系（玻璃、晶体），分析网络形成体 (former) 与配体之间的连接拓扑。

```bash
# P-O 体系（截断半径 2.3 Å）
mol-network -i traj.dump --P-O=2.3

# 多元素体系，输出 xlsx
mol-network -i traj.dump --P-O=2.3 --Si-O=1.8 --Al-O=2.0 --format xlsx -o result

# 参数格式：--Former-Ligand=cutoff_in_Angstrom
# Former 首字母大写（与 clap 普通参数区分）
```

输出内容：
- **CN 分布**：每种 (former, ligand) 对的配位数分布 + 平均值
- **配体分类**：每个配体原子的类型（FO / NBO(X) / BO(X-Y) / OBO(X-Y-Z)）
- **Qn 物种**：Q^n(mX) 分布（n = 桥接配体数，m = 异元素桥接数）

CSV 输出三个文件（`_cn.csv`、`_ligand.csv`、`_qn.csv`）；xlsx 输出三个 Sheet。

## 内部单位

| 物理量 | 单位 |
|---|---|
| 长度 | Å |
| 能量 | eV |
| 力 | eV/Å |
| 应力 | eV/Å³ |
| 时间 | fs |
| 质量 | amu |

## 核心数据类型

```rust
// 始终以 Trajectory 作为顶层类型，单帧文件也是 Trajectory { frames: vec![frame] }
pub struct Trajectory { pub frames: Vec<Frame>, pub metadata: TrajectoryMetadata }

pub struct Frame {
    pub atoms: Vec<Atom>,
    pub cell: Option<Cell>,           // None = 非周期性
    pub pbc: [bool; 3],               // 各轴是否启用 PBC
    pub energy: Option<f64>,
    pub forces: Option<Vec<Vector3<f64>>>,
    pub velocities: Option<Vec<Vector3<f64>>>,
    // ...
}

pub struct Atom {
    pub element: String,
    pub position: Vector3<f64>,       // Å，直角坐标
    pub label: Option<String>,
    pub mass: Option<f64>,
    pub magmom: Option<f64>,
    pub charge: Option<f64>,
}
```

## 常用开发命令

```bash
cargo build                          # 编译工作区
cargo build --release                # 发布版
cargo test                           # 全部测试（当前 146 个）
cargo test --package molflow-io      # 单个 crate 测试
cargo fmt && cargo clippy            # 格式化 + 静态检查
cargo run --bin mol-info -- -i examples/water.xyz
```

## 测试文件

`tests/`（water.xyz, water.pdb, water.cif）、`examples/`（LAMMPS dump、CIF、LMP data、CP2K inp）
