# ChemTools 快速入门指南

## 新手学习路线

### 第一步：理解项目结构
```
molflow/
├── molflow-core/      ← 从这里开始！定义了 Atom, Molecule, Trajectory
├── molflow-io/        ← 文件读写：read_xyz(), write_pdb()
├── molflow-analysis/  ← 分析功能：计算距离、RMSD 等
├── molflow-workflow/  ← 生成输入文件：Gaussian, GROMACS
├── molflow-cli/       ← 命令行工具：你将主要使用这个
└── molflow-python/    ← Python 绑定：方便与 Python 集成
```

### 第二步：编译项目
```bash
# 在 molflow 根目录
cargo build --release

# 如果成功，你会看到编译输出
# 可执行文件在: target/release/molflow
```

### 第三步：运行第一个命令
```bash
# 查看帮助
cargo run --bin molflow -- --help

# 显示水分子信息
cargo run --bin molflow -- info -i examples/water.xyz

# 输出应该显示:
# 分子信息:
#   原子数: 3
#   总质量: 18.02 amu
#   质心: (0.000, 0.000, 0.065)
```

### 第四步：理解代码流程

当你运行 `molflow info -i water.xyz` 时，发生了什么？

1. **CLI 解析** (`molflow-cli/src/main.rs`)
   - clap 解析命令行参数
   - 调用 `show_info()` 函数

2. **读取文件** (`molflow-io/src/readers/xyz.rs`)
   - `read_xyz()` 打开文件
   - 解析原子坐标
   - 创建 `Molecule` 对象

3. **使用核心** (`molflow-core/src/molecule.rs`)
   - 调用 `mol.atom_count()`
   - 调用 `mol.total_mass()`
   - 调用 `mol.center_of_mass()`

4. **输出结果**
   - 打印到终端

## 常用操作示例

### 1. 文件格式转换
```bash
# XYZ -> PDB
cargo run --bin molflow -- convert -i examples/water.xyz -o water.pdb

# 查看生成的 PDB 文件
cat water.pdb
```

### 2. 创建 Gaussian 输入文件
```bash
cargo run --bin molflow -- job \
    -i examples/water.xyz \
    -s gaussian \
    -m B3LYP \
    -o water.gjf

# 查看生成的任务文件
cat water.gjf
```

### 3. 在 Rust 代码中使用
创建文件 `examples/my_first_program.rs`:

```rust
use molflow_io::read_xyz;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // 读取分子
    let mol = read_xyz("examples/water.xyz")?;
    
    // 打印信息
    println!("分子有 {} 个原子", mol.atom_count());
    
    // 遍历所有原子
    for (i, atom) in mol.atoms.iter().enumerate() {
        println!("原子 {}: {} 在 ({:.3}, {:.3}, {:.3})",
            i + 1,
            atom.element,
            atom.position.x,
            atom.position.y,
            atom.position.z
        );
    }
    
    Ok(())
}
```

运行：
```bash
# 在项目根目录创建示例程序
cargo new --bin my_example
cd my_example
# 在 Cargo.toml 中添加依赖
# [dependencies]
# molflow-io = { path = "../molflow-io" }
# molflow-core = { path = "../molflow-core" }
# nalgebra = "0.32"
# anyhow = "1.0"

cargo run
```

## 如何添加新功能

### 示例：添加新的文件格式 (MOL2)

1. **在 molflow-io 中添加读取器**
   ```rust
   // molflow-io/src/readers/mol2.rs
   pub fn read_mol2(path: &str) -> Result<Molecule> {
       // 实现 MOL2 解析逻辑
       todo!()
   }
   ```

2. **在 mod.rs 中导出**
   ```rust
   // molflow-io/src/readers/mod.rs
   pub mod mol2;
   pub use mol2::read_mol2;
   ```

3. **在 CLI 中添加支持**
   ```rust
   // 在 convert_file() 函数中添加格式判断
   match input.split('.').last() {
       Some("mol2") => read_mol2(input)?,
       Some("xyz") => read_xyz(input)?,
       // ...
   }
   ```

## 学习建议

### 对于 Rust 新手：
1. 先学习基本的 Rust 语法
2. 理解所有权（Ownership）概念
3. 从阅读 `molflow-core` 开始，它最简单
4. 尝试修改现有代码，看看会发生什么
5. 编译器是你的朋友！仔细读错误信息

### 对于计算化学专业人员：
1. 关注 `molflow-io` 和 `molflow-workflow`
2. 这些模块直接对应你熟悉的文件格式和软件
3. 从添加你常用的文件格式开始
4. 逐步添加你的工作流程

### 推荐的学习顺序：
1. 运行现有的命令行工具，熟悉功能
2. 阅读 `molflow-core` 的代码，理解数据结构
3. 阅读 `molflow-io` 的 XYZ 读取器，理解文件解析
4. 尝试修改 CLI 工具，添加新命令
5. 实现一个新的文件格式支持
6. 添加自己需要的分析功能

## 调试技巧

```bash
# 查看详细编译信息
cargo build --verbose

# 运行测试
cargo test

# 查看某个模块的测试
cargo test --package molflow-core

# 格式化代码
cargo fmt

# 检查代码（比 build 快）
cargo check

# 查看依赖树
cargo tree
```

## 遇到问题？

1. **编译错误**：仔细阅读错误信息，Rust 编译器会告诉你怎么修复
2. **找不到函数**：检查 `use` 语句是否正确导入
3. **类型不匹配**：查看函数签名，确保类型正确
4. **文件找不到**：确保路径正确，使用相对路径或绝对路径

## 下一步

- 尝试实现一个新的分析功能（如计算原子间距离）
- 添加对你常用文件格式的支持
- 将项目集成到你的工作流程中
- 为项目贡献代码！
