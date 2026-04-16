# ChemTools - 计算化学工具集

一个用 Rust 编写的模块化计算化学工具包，支持分子结构处理、轨迹分析和计算任务创建。

## 项目结构

```
molflow/
├── molflow-core/          # 核心数据结构
├── molflow-io/            # 文件读写
├── molflow-analysis/      # 后处理分析
├── molflow-workflow/      # 任务创建
├── molflow-cli/           # 命令行工具
└── molflow-python/        # Python 绑定
```

## 模块说明

### 1. molflow-core - 核心模块
**作用**: 定义基础数据结构
- `Atom`: 原子结构（位置、元素、质量、电荷等）
- `Molecule`: 分子结构（原子列表、键连接）
- `Trajectory`: 分子动力学轨迹（多帧分子结构）
- `units`: 单位转换工具（长度、能量、温度）
- `error`: 统一错误处理

**使用场景**: 所有其他模块的基础，提供统一的数据表示

### 2. molflow-io - 文件读写模块
**作用**: 支持多种文件格式的读取和写入
- **读取器**: XYZ, PDB, GRO（可扩展）
- **写入器**: XYZ, PDB, GRO（可扩展）

**使用场景**: 
- 读取实验或计算得到的结构文件
- 导出结果供其他软件使用
- 格式转换

### 3. molflow-analysis - 分析模块
**作用**: 提供后处理分析功能
- **几何分析**: 键长、键角、二面角、回转半径
- **轨迹分析**: RMSD, MSD, 质心轨迹
- **性质计算**: 偶极矩、能量等

**使用场景**:
- 分析分子动力学轨迹
- 提取几何参数
- 计算分子性质

### 4. molflow-workflow - 任务创建模块
**作用**: 生成各种量子化学软件的输入文件
- **支持软件**: Gaussian, ORCA, GROMACS（可扩展）
- **任务类型**: 几何优化、单点能、频率计算、MD模拟
- **模板管理**: 常用计算任务模板

**使用场景**:
- 批量生成计算任务
- 标准化工作流程
- 快速构建输入文件

### 5. molflow-cli - 命令行工具
**作用**: 提供便捷的命令行接口
- 文件格式转换
- 快速分析
- 任务创建
- 信息查看

**使用场景**: 
- 日常快速操作
- Shell 脚本集成
- 自动化工作流

### 6. molflow-python - Python 绑定
**作用**: 在 Python 中使用 Rust 核心功能
- 高性能计算核心
- Python 易用接口
- 与其他 Python 库集成

**使用场景**:
- Jupyter Notebook 交互分析
- 与 NumPy/Pandas 集成
- Python 脚本自动化

## 快速开始

### 安装 Rust (如果还没有)
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

### 编译项目
```bash
cd molflow
cargo build --release
```

### 使用命令行工具
```bash
# 文件格式转换
cargo run --bin molflow -- convert -i input.xyz -o output.pdb

# 显示分子信息
cargo run --bin molflow -- info -i molecule.xyz

# 创建 Gaussian 任务
cargo run --bin molflow -- job -i input.xyz -s gaussian -m B3LYP -o job.gjf
```

### 作为 Rust 库使用
```rust
use molflow_core::{Molecule, Atom};
use molflow_io::{read_xyz, write_pdb};
use nalgebra::Vector3;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // 读取文件
    let mol = read_xyz("input.xyz")?;
    
    // 计算质心
    let com = mol.center_of_mass();
    println!("质心: {:?}", com);
    
    // 导出为 PDB
    write_pdb(&mol, "output.pdb")?;
    
    Ok(())
}
```

### Python 使用 (可选，需要先安装 maturin)

⚠️ **新手建议**: 先跳过 Python 部分，专注于 Rust 核心功能。

详细说明请看 `PYTHON_BINDING.md`

```bash
cd molflow-python
pip install maturin
maturin develop
```

```python
import molflow

# 读取文件
info = molflow.read_xyz("molecule.xyz")
print(info)

# 计算质心
com = molflow.center_of_mass("molecule.xyz")
print(f"质心: {com}")
```

## 扩展指南

### 添加新的文件格式
1. 在 `molflow-io/src/readers/` 添加新的读取器
2. 在 `molflow-io/src/writers/` 添加新的写入器
3. 在 `mod.rs` 中导出

### 添加新的分析功能
1. 在 `molflow-analysis/src/` 添加新模块
2. 实现分析函数
3. 在 CLI 和 Python 绑定中添加接口

### 添加新的计算软件支持
1. 在 `molflow-workflow/src/job_builder.rs` 添加新的构建器
2. 在 `templates.rs` 添加模板
3. 在 CLI 中添加命令

## 设计优势

### 为什么选择 Workspace 结构？
1. **模块化**: 每个功能独立开发、测试
2. **可复用**: 其他项目可以只依赖需要的模块
3. **清晰职责**: core → io → analysis → workflow 层次分明
4. **易维护**: 修改一个模块不影响其他模块
5. **并行开发**: 团队可以同时开发不同模块

### 为什么使用 Rust？
1. **性能**: 接近 C/C++ 的性能
2. **安全**: 编译时内存安全保证
3. **Cargo**: 无需学习 Make/CMake
4. **生态**: 丰富的科学计算库
5. **Python 互操作**: 通过 PyO3 轻松集成

## 开发路线图

- [x] 基础框架搭建
- [ ] 完善文件格式支持 (MOL2, SDF, LAMMPS 等)
- [ ] 实现完整的轨迹分析功能
- [ ] 支持更多量子化学软件 (ORCA, Q-Chem, CP2K 等)
- [ ] 添加周期性边界条件支持
- [ ] 实现力场参数处理
- [ ] 添加可视化功能
- [ ] 性能优化和并行化

## 依赖说明

- `nalgebra`: 线性代数（向量、矩阵）
- `ndarray`: 多维数组
- `rayon`: 并行计算
- `serde`: 序列化/反序列化
- `thiserror`: 错误处理
- `anyhow`: 错误传播
- `clap`: 命令行参数解析
- `pyo3`: Python 绑定

## 许可证

MIT License

## 贡献指南

欢迎提交 Issue 和 Pull Request！

对于新手：
1. 从简单的功能开始（如添加新的文件格式）
2. 保持代码简洁、注释清晰
3. 添加单元测试
4. 更新文档
