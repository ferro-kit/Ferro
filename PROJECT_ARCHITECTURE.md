# ChemTools 项目架构完整说明

## 🎯 项目概述

这是一个为计算化学领域设计的 Rust 工具包，采用 **Workspace 多包架构**，将功能模块化，便于长期维护和扩展。

## 📦 模块架构图

```
┌─────────────────────────────────────────────────────────┐
│                    molflow-cli                         │
│                  (命令行工具入口)                          │
│  cargo run --bin molflow -- convert/analyze/job        │
└───────────┬────────────┬────────────┬───────────────────┘
            │            │            │
            ▼            ▼            ▼
┌───────────────┐  ┌───────────┐  ┌──────────────────┐
│molflow-io   │  │molflow- │  │ molflow-       │
│(文件读写)      │  │analysis   │  │ workflow         │
│- read_xyz()   │  │(后处理)    │  │ (任务创建)        │
│- write_pdb()  │  │- RMSD     │  │ - Gaussian       │
│- read_pdb()   │  │- 几何计算  │  │ - GROMACS        │
└───────┬───────┘  └─────┬─────┘  └────────┬─────────┘
        │                │                  │
        │                │                  │
        └────────────────┼──────────────────┘
                         ▼
              ┌──────────────────┐
              │ molflow-core   │
              │  (核心数据结构)    │
              │  - Atom          │
              │  - Molecule      │
              │  - Trajectory    │
              └──────────────────┘
                         │
                         │ (可选)
                         ▼
              ┌──────────────────┐
              │ molflow-python │
              │   (Python绑定)    │
              │   import molflow│
              └──────────────────┘
```

## 🔧 每个模块的作用

### 1. molflow-core (核心基础)
**职责**: 定义所有基础数据类型
- `Atom`: 单个原子（元素、坐标、质量、电荷）
- `Molecule`: 分子（原子列表、化学键）
- `Trajectory`: 轨迹（时间序列的分子结构）
- `units`: 单位转换（Å↔nm, Hartree↔kcal/mol）
- `error`: 统一错误类型

**何时修改**: 需要新的数据结构时（如添加周期性盒子）

### 2. molflow-io (输入输出)
**职责**: 文件格式的读取和写入
- **readers/**: XYZ, PDB, GRO 读取器
- **writers/**: XYZ, PDB, GRO 写入器

**何时修改**: 
- 添加新格式（MOL2, SDF, LAMMPS data）
- 优化解析性能
- 支持压缩文件

### 3. molflow-analysis (数据分析)
**职责**: 后处理计算和分析
- **geometry**: 键长、键角、二面角
- **trajectory_analysis**: RMSD, MSD, 回转半径
- **properties**: 偶极矩、能量等分子性质

**何时修改**:
- 添加新分析方法（如径向分布函数）
- 实现统计分析
- 添加可视化功能

### 4. molflow-workflow (任务管理)
**职责**: 生成计算软件的输入文件
- **job_builder**: Gaussian, ORCA, GROMACS 任务构建器
- **templates**: 常用任务模板

**何时修改**:
- 支持新软件（CP2K, Quantum ESPRESSO）
- 添加新任务类型（QM/MM, 过渡态搜索）
- 实现批量任务生成

### 5. molflow-cli (命令行工具)
**职责**: 提供用户友好的命令行界面
- `convert`: 文件格式转换
- `analyze`: 轨迹分析
- `job`: 创建计算任务
- `info`: 显示分子信息

**何时修改**:
- 添加新命令
- 改进用户体验
- 添加交互式模式

### 6. molflow-python (Python 集成)
**职责**: 允许 Python 调用 Rust 核心功能
- 提供 Python API
- 保持高性能核心
- 与 NumPy/Pandas 互操作

**何时修改**:
- 需要在 Jupyter 中使用
- 与 Python 生态集成
- 创建 Web 应用

## 🚀 使用流程

### 场景 1: 日常文件转换
```bash
# 用户操作
molflow convert -i structure.xyz -o structure.pdb

# 内部流程
CLI (main.rs) → molflow-io (read_xyz) → 
molflow-core (Molecule) → molflow-io (write_pdb)
```

### 场景 2: 创建计算任务
```bash
# 用户操作
molflow job -i mol.xyz -s gaussian -m B3LYP

# 内部流程
CLI → molflow-io (read_xyz) → 
molflow-workflow (GaussianJobBuilder) → 生成 .gjf 文件
```

### 场景 3: Python 脚本分析
```python
import molflow
mol = molflow.read_xyz("traj.xyz")
com = molflow.center_of_mass("traj.xyz")

# 内部流程
Python → PyO3 → molflow-io → molflow-core → 返回结果
```

## 📝 开发工作流

### 添加新功能的标准流程

1. **明确需求**: 确定功能属于哪个模块
2. **修改 core**: 如需新数据结构，先修改 core
3. **实现功能**: 在对应模块添加代码
4. **添加测试**: 写单元测试确保正确性
5. **更新 CLI**: 如需暴露给用户，更新 CLI
6. **更新文档**: 更新 README 和注释

### 示例: 添加 MOL2 格式支持

```
1. molflow-io/src/readers/mol2.rs  ← 创建读取器
2. molflow-io/src/writers/mol2.rs  ← 创建写入器
3. molflow-io/src/readers/mod.rs   ← 导出模块
4. molflow-cli/src/main.rs         ← 添加格式判断
5. 测试: cargo test --package molflow-io
6. 更新 README.md
```

## 🎓 新手学习路径

### 第 1 周: 熟悉项目
- 阅读 README.md 和 QUICKSTART.md
- 编译项目: `cargo build`
- 运行示例: `cargo run --bin molflow -- info -i examples/water.xyz`
- 理解数据流: CLI → IO → Core

### 第 2 周: 阅读代码
- 从 `molflow-core` 开始，理解 Atom 和 Molecule
- 看 `molflow-io` 的 XYZ 读取器
- 理解 `molflow-cli` 的命令分发

### 第 3 周: 简单修改
- 在 Atom 中添加一个新字段
- 修改 XYZ 读取器支持注释中的额外信息
- 添加一个新的 CLI 命令

### 第 4 周: 独立功能
- 实现一个新的文件格式支持
- 添加一个新的分析功能
- 为新功能编写测试

## 🔍 依赖关系

```
molflow-cli
    ├── molflow-core
    ├── molflow-io
    │   └── molflow-core
    ├── molflow-analysis
    │   └── molflow-core
    └── molflow-workflow
        └── molflow-core

molflow-python
    ├── molflow-core
    ├── molflow-io
    └── molflow-analysis
```

**设计原则**:
- core 不依赖任何其他模块（最底层）
- 其他模块只依赖 core
- CLI 和 Python 绑定依赖所有功能模块

## 🛠️ 常用命令

```bash
# 编译整个项目
cargo build --release

# 只编译某个模块
cargo build --package molflow-core

# 运行测试
cargo test

# 运行特定模块的测试
cargo test --package molflow-io

# 格式化代码
cargo fmt

# 代码检查（不编译）
cargo check

# 查看依赖树
cargo tree

# 更新依赖
cargo update

# 清理编译产物
cargo clean
```

## 📊 项目统计

- 总模块数: 6
- 核心模块: 1 (molflow-core)
- 功能模块: 3 (io, analysis, workflow)
- 接口模块: 2 (cli, python)
- 编程语言: Rust + Python (绑定)
- 依赖管理: Cargo (Workspace)

## 🎯 为什么选择这个架构？

1. **模块化**: 每个模块职责单一，易于理解
2. **可扩展**: 添加新功能不影响现有代码
3. **可复用**: 其他项目可以只使用需要的模块
4. **易测试**: 每个模块独立测试
5. **并行开发**: 团队可以同时开发不同模块
6. **性能优化**: 只重编译修改的模块

## 📌 注意事项

1. **不要跨层依赖**: analysis 不应该直接依赖 io
2. **保持 core 纯净**: core 只包含数据结构，不含 I/O 逻辑
3. **错误处理**: 库使用 Result<T, ChemError>，CLI 使用 anyhow
4. **文档注释**: 所有公开函数都应有 /// 注释
5. **单元测试**: 每个模块都应有测试

## 🚀 未来扩展方向

1. **更多文件格式**: MOL2, SDF, LAMMPS, AMBER
2. **高级分析**: 自由能计算、溶剂化分析
3. **可视化**: 集成 3D 可视化
4. **并行优化**: 大规模轨迹并行处理
5. **机器学习**: 集成分子描述符和预测模型
6. **Web 界面**: 基于 Python 绑定的 Web 应用

## 📖 参考资源

- Rust 官方文档: https://doc.rust-lang.org/book/
- Cargo Workspace: https://doc.rust-lang.org/cargo/reference/workspaces.html
- PyO3 文档: https://pyo3.rs/
- nalgebra 文档: https://nalgebra.org/

祝你开发愉快！🎉
