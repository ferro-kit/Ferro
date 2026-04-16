# Python 绑定使用指南

## ⚠️ 新手建议

如果你是 Rust 新手，**建议先跳过 Python 绑定**，专注于：
1. 学习 Rust 基础
2. 熟悉项目结构
3. 使用命令行工具
4. 添加新功能

等熟悉了 Rust 和项目后，再回来配置 Python 绑定。

## 为什么 `cargo build` 会失败？

Python 扩展模块需要链接到 Python 解释器，不能用普通的 `cargo build`。
必须使用专门的构建工具：`maturin`

## 如何启用 Python 绑定

### 第 1 步：确保 Python 已安装

```bash
# 检查 Python 版本（需要 3.8+）
python3 --version

# 或
python --version
```

### 第 2 步：安装 maturin

```bash
pip install maturin
# 或
pip3 install maturin
```

### 第 3 步：取消注释 Cargo.toml

编辑根目录的 `Cargo.toml`，取消注释：

```toml
[workspace]
members = [
    "molflow-core",
    "molflow-io",
    "molflow-analysis",
    "molflow-workflow",
    "molflow-cli",
    "molflow-python",  # 去掉前面的 #
]
```

### 第 4 步：构建 Python 模块

```bash
cd molflow-python

# 开发模式（安装到当前 Python 环境）
maturin develop

# 或构建 wheel 包
maturin build --release
```

### 第 5 步：使用

```python
import molflow

# 读取文件
info = molflow.read_xyz("../examples/water.xyz")
print(info)

# 计算质心
com = molflow.center_of_mass("../examples/water.xyz")
print(f"质心: {com}")
```

## 常见问题

### 问题 1: "Python.h not found"

**原因**: 缺少 Python 开发头文件

**解决方案**:
```bash
# macOS
brew install python3

# Ubuntu/Debian
sudo apt-get install python3-dev

# CentOS/RHEL
sudo yum install python3-devel
```

### 问题 2: 链接错误（symbol not found）

**原因**: 
- Python 版本不匹配
- 需要使用 maturin 而不是 cargo

**解决方案**:
```bash
# 清理旧的构建
cargo clean

# 使用 maturin
cd molflow-python
maturin develop
```

### 问题 3: 多个 Python 版本冲突

**解决方案**:
```bash
# 指定 Python 解释器
maturin develop --interpreter python3.11

# 或使用虚拟环境
python3 -m venv venv
source venv/bin/activate
maturin develop
```

## 不使用 Python 绑定也能完成工作

记住：Python 绑定是**可选的**。你可以：

### 方案 1: 只用 Rust CLI
```bash
# 格式转换
molflow convert -i input.xyz -o output.pdb

# 分析
molflow analyze -t traj.xyz -p rmsd

# 创建任务
molflow job -i mol.xyz -s gaussian
```

### 方案 2: 用 Rust 写脚本
```rust
// my_script.rs
use molflow_io::{read_xyz, write_pdb};

fn main() {
    let mol = read_xyz("input.xyz").unwrap();
    write_pdb(&mol, "output.pdb").unwrap();
    println!("完成!");
}
```

### 方案 3: 编译后用 shell 调用
```bash
#!/bin/bash
# process_all.sh
for file in *.xyz; do
    molflow convert -i "$file" -o "${file%.xyz}.pdb"
done
```

## 什么时候需要 Python 绑定？

只有在以下情况下才真正需要：
- ✅ 需要在 Jupyter Notebook 中交互分析
- ✅ 需要与 NumPy/Pandas/Matplotlib 集成
- ✅ 需要在 Python web 应用中使用
- ✅ 团队主要使用 Python

如果你的主要工作是：
- ❌ 批量处理文件 → 用 CLI
- ❌ 自动化工作流 → 用 CLI + shell
- ❌ 学习 Rust → 直接写 Rust

## 推荐学习路径

1. **第 1-2 周**: 熟悉 Rust 核心功能
   - 使用 CLI 工具
   - 阅读代码
   - 添加简单功能

2. **第 3-4 周**: 扩展功能
   - 添加新文件格式
   - 实现分析功能
   - 优化性能

3. **第 5+ 周**: （可选）Python 集成
   - 配置 maturin
   - 构建 Python 模块
   - 在 Jupyter 中使用

## 总结

- 🎯 **新手**: 先注释掉 Python 模块，专注 Rust
- 🔧 **中级**: 熟悉项目后，使用 maturin 构建
- 🚀 **高级**: 根据需求决定是否需要 Python 绑定

记住：**Python 绑定不是必需品，是锦上添花！**
