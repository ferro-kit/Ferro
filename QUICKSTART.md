# molflow 快速入门

## 编译

```bash
cd molflow
cargo build --release
# 可执行文件在 target/release/mol-*
```

或直接用 `cargo run --bin <name> -- <args>` 运行（无需先 build）。

---

## 常用操作速查

### 格式转换
```bash
mol-convert -i structure.xyz -o structure.pdb
mol-convert -i traj.dump    -o traj.extxyz
mol-convert -i NaCl.cif     -o POSCAR
```

### 查看结构信息
```bash
mol-info -i water.xyz
mol-info -i NaCl.cif
mol-info -i traj.dump        # 显示帧数、原子数、元素组成、晶胞参数
```

### 生成 Gaussian 输入文件
```bash
mol-job -i water.xyz -s gaussian -m B3LYP -b "6-31G*" -o job.gjf
```

---

## 轨迹分析

所有分析命令：无 `-i` 时打印该模式的详细帮助（用途 / 参数 / 示例）。

### 径向分布函数
```bash
mol-traj -m gr -i traj.dump
mol-traj -m gr -i traj.dump --r-max 8.0 --r-cut 3.2 --last-n 500
# 输出：gr.dat + gr_cn.dat
```

### 结构因子
```bash
mol-traj -m sq -i traj.dump
mol-traj -m sq -i traj.dump --weighting xrd --q-max 20.0
```

### 均方位移
```bash
mol-traj -m msd -i traj.dump --dt 2.0
mol-traj -m msd -i traj.dump --elements Li --dt 1.0 --last-n 2000
```

### 键角分布
```bash
mol-traj -m angle -i traj.xyz --r-cut-ab 2.0
```

---

## 相关函数

```bash
mol-corr -m vacf    -i traj.dump --dt 2.0        # 速度自相关
mol-corr -m rotcorr -i traj.xyz --center O --neighbor H   # 转动相关
mol-corr -m vanhove -i traj.dump --tau 100        # Van Hove 自相关
```

---

## 空间分布（cube 格式）

```bash
mol-cube -m density  -i traj.dump --nx 100 --ny 100 --nz 100
mol-cube -m velocity -i traj.dump --elements Li
mol-cube -m force    -i traj.dump -o force.cube
# 输出为 Gaussian cube，用 VESTA / VMD 打开
```

---

## 玻璃网络结构分析

分析配位数、配体类型（FO/NBO/BO/OBO）和 Qn 物种分布。

```bash
# 磷酸盐体系
mol-network -i traj.dump --P-O=2.3

# 多元素体系
mol-network -i traj.dump --P-O=2.3 --Si-O=1.8 --Al-O=2.0

# 含 F 的体系
mol-network -i traj.dump --P-O=2.3 --P-F=2.1

# 输出 xlsx（三个 Sheet：CN / Ligand / Qn）
mol-network -i traj.dump --P-O=2.3 --format xlsx -o result

# 只分析最后 500 帧
mol-network -i traj.dump --P-O=2.3 --last-n 500
```

参数格式：`--Former-Ligand=cutoff_Å`
- Former 首字母大写（与 `--last-n` 等普通参数区分）
- 无 `-i` 时打印使用帮助

---

## 并行与性能

所有分析命令支持 `--ncore N` 指定线程数（默认使用全部核心）。

```bash
mol-traj -m gr -i large_traj.dump --ncore 8
mol-cube -m density -i traj.dump --ncore 16
```

---

## 测试

```bash
cargo test                          # 全部测试（约 146 个）
cargo test --package molflow-io     # 单 crate
cargo test --package molflow-analysis network   # 指定模块
```

---

## 调试编译

```bash
cargo check                  # 快速语法检查（不生成二进制）
cargo build                  # 调试版（更快编译，较慢运行）
cargo build --release        # 发布版（慢编译，快运行）
cargo fmt && cargo clippy    # 格式化 + lint
```
