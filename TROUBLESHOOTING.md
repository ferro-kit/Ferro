# ChemTools 常见问题解决

## 编译错误

### 错误 1: nalgebra Vector3 不支持 serde 序列化

**错误信息**:
```
error[E0277]: the trait bound `Matrix<f64, Const<3>, Const<1>, ArrayStorage<f64, 3, 1>>: serde::Serialize` is not satisfied
```

**原因**: 
nalgebra 的 `Vector3<f64>` 类型默认不支持 serde 序列化。需要启用 nalgebra 的 `serde-serialize` 特性。

**解决方案**:
在 `Cargo.toml` 中修改 nalgebra 依赖：

```toml
# ❌ 错误写法
nalgebra.workspace = true

# ✅ 正确写法
nalgebra = { workspace = true, features = ["serde-serialize"] }
```

**已修复的文件**:
- `molflow-core/Cargo.toml`
- `molflow-analysis/Cargo.toml`
- `molflow-io/Cargo.toml`

---

### 错误 2: 缺少依赖导致编译失败

**错误信息**:
```
error[E0432]: unresolved import `nalgebra`
 --> molflow-io/src/readers/xyz.rs:6:5
  |
6 | use nalgebra::Vector3;
  |     ^^^^^^^^ use of unresolved module or unlinked crate `nalgebra`
```

**原因**: 
在代码中使用了某个 crate，但是在 `Cargo.toml` 中没有声明依赖。

**解决方案**:
在对应模块的 `Cargo.toml` 中添加依赖：

```toml
[dependencies]
molflow-core = { path = "../molflow-core" }
nalgebra = { workspace = true, features = ["serde-serialize"] }  # 添加这一行
```

**已修复**: `molflow-io` 现在已经正确声明了 `nalgebra` 依赖。

---

### 错误 3: Python 绑定链接失败

**错误信息**:
```
error: linking with `cc` failed: exit status: 1
...
Undefined symbols for architecture arm64:
  "_PyBytes_AsString", referenced from:
  ...
ld: symbol(s) not found for architecture arm64
```

**原因**: 
Python 扩展模块不能用 `cargo build` 构建，必须使用 `maturin` 工具。

**解决方案 1: 跳过 Python 模块（推荐新手）**

编辑根目录 `Cargo.toml`：
```toml
[workspace]
members = [
    "molflow-core",
    "molflow-io",
    "molflow-analysis",
    "molflow-workflow",
    "molflow-cli",
    # "molflow-python",  # 添加 # 注释掉这一行
]
```

然后重新编译：
```bash
cargo build --release
```

**解决方案 2: 正确构建 Python 模块**

```bash
# 安装 maturin
pip install maturin

# 构建 Python 模块
cd molflow-python
maturin develop
```

**详细说明**: 查看 `PYTHON_BINDING.md` 文档

**建议**: 作为新手，先专注于 Rust 部分，等熟悉后再研究 Python 绑定。

---

## 其他常见问题

### 问题 2: 找不到 cargo 命令

**解决方案**: 确保已安装 Rust
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env
```

### 问题 3: 编译速度慢

**解决方案**: 使用增量编译和并行构建
```bash
# 只检查代码，不生成可执行文件（快）
cargo check

# 只编译修改的模块
cargo build --package molflow-core

# 使用更多 CPU 核心
cargo build -j 8
```

### 问题 4: 依赖下载失败

**解决方案**: 
```bash
# 使用国内镜像（中国用户）
# 在 ~/.cargo/config.toml 添加：
[source.crates-io]
replace-with = 'ustc'

[source.ustc]
registry = "https://mirrors.ustc.edu.cn/crates.io-index"
```

### 问题 5: 运行测试失败

**解决方案**:
```bash
# 显示详细测试输出
cargo test -- --nocapture

# 运行特定测试
cargo test test_atom_creation

# 只测试某个包
cargo test --package molflow-core
```

---

## Rust 学习中的常见困惑

### 所有权 (Ownership)

```rust
// ❌ 错误：值已被移动
let mol = Molecule::new();
process(mol);  // mol 被移动
println!("{:?}", mol);  // 编译错误！

// ✅ 方法 1: 使用引用
process(&mol);  // 借用，不移动
println!("{:?}", mol);  // 正确

// ✅ 方法 2: 克隆
process(mol.clone());  // 克隆一份
println!("{:?}", mol);  // 正确
```

### Option 类型

```rust
// ❌ 可能 panic
let atom = mol.atoms[100];  // 如果不存在会 panic

// ✅ 安全的方式
if let Some(atom) = mol.atoms.get(100) {
    // 使用 atom
}

// ✅ 或者使用 Option 的方法
mol.atoms.get(100).map(|atom| {
    // 使用 atom
});
```

### Result 错误处理

```rust
// ❌ 不处理错误
let mol = read_xyz("file.xyz").unwrap();  // 可能 panic

// ✅ 使用 ? 操作符
fn process() -> Result<(), Box<dyn Error>> {
    let mol = read_xyz("file.xyz")?;  // 自动传播错误
    Ok(())
}

// ✅ 或者 match
match read_xyz("file.xyz") {
    Ok(mol) => println!("成功"),
    Err(e) => println!("错误: {}", e),
}
```

---

## 特性标志 (Feature Flags) 说明

### 什么是特性标志？

Rust 的 crate 可以有可选功能，通过特性标志启用。这样可以：
- 减小编译时间和二进制大小
- 按需启用功能
- 避免不必要的依赖

### 常见特性

```toml
# nalgebra 的常用特性
nalgebra = { version = "0.32", features = [
    "serde-serialize",  # serde 支持
    "rayon",            # 并行计算
    "mint",             # 与其他数学库互操作
] }

# serde 的常用特性
serde = { version = "1.0", features = [
    "derive",  # 使用 #[derive(Serialize, Deserialize)]
] }

# 查看 crate 支持的特性
cargo tree -e features
```

---

## 调试技巧

### 1. 使用 dbg! 宏

```rust
let mol = read_xyz("file.xyz")?;
dbg!(&mol.atom_count());  // 打印调试信息
```

### 2. 添加打印语句

```rust
println!("调试: mol 有 {} 个原子", mol.atom_count());
```

### 3. 使用 Rust Analyzer (VS Code)

- 安装 rust-analyzer 插件
- 鼠标悬停查看类型
- F12 跳转到定义
- Ctrl+Space 自动补全

### 4. 查看编译错误的详细信息

```bash
cargo build --verbose
rustc --explain E0277  # 查看错误码说明
```

---

## 性能优化建议

### 1. Release 模式

```bash
# Debug 模式（慢，但有调试信息）
cargo build

# Release 模式（快，优化后）
cargo build --release
```

### 2. 使用并行计算

```rust
use rayon::prelude::*;

// 串行
for frame in &traj.frames {
    process(frame);
}

// 并行
traj.frames.par_iter().for_each(|frame| {
    process(frame);
});
```

### 3. 避免不必要的克隆

```rust
// ❌ 慢：每次都克隆
for atom in &mol.atoms {
    let a = atom.clone();  // 不必要
    process(&a);
}

// ✅ 快：直接使用引用
for atom in &mol.atoms {
    process(atom);
}
```

---

## 获取帮助

1. **Rust 官方文档**: https://doc.rust-lang.org/book/
2. **Rust by Example**: https://doc.rust-lang.org/rust-by-example/
3. **Rust 社区**: https://users.rust-lang.org/
4. **Stack Overflow**: 搜索 [rust] + 你的问题
5. **查看 crate 文档**: https://docs.rs/

---

## 项目特定问题

### 如何添加新的文件格式？

参考 `molflow-io/src/readers/xyz.rs`，按照相同模式创建。

### 如何添加新的分析功能？

在 `molflow-analysis` 中添加新函数，然后在 CLI 中添加命令。

### 如何为 Python 添加新函数？

在 `molflow-python/src/lib.rs` 中添加 `#[pyfunction]`。

---

**提示**: 编译器是你的朋友！Rust 的错误信息非常详细，仔细阅读通常能找到解决方案。
