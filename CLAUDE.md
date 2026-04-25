# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Test Commands

```bash
cargo build                                                    # build entire workspace
cargo build --release
cargo build --package nexflux-core                             # single crate
cargo test                                                     # all tests
cargo test --package nexflux-io                                # single crate
cargo test --package nexflux-core test_basic_molecule          # single test
cargo fmt
cargo clippy
cargo check

# CLI (dev)
cargo run --bin nex-convert -- -i input.xyz -o output.pdb
cargo run --bin nex-info    -- -i input.xyz
cargo run --bin nex-job     -- -i input.xyz -s gaussian -m B3LYP -o job.gjf
cargo run --bin nex-traj    -- -m gr  -i traj.dump
cargo run --bin nex-cube    -- -m sdf -i traj.dump --qn 3 --modifier Zn

# Python bindings (requires maturin; nexflux-python not yet in workspace)
cd nexflux-python && maturin develop
```

Test fixtures: `tests/` (water.xyz, water.pdb, water.cif)

---

## Architecture

Cargo workspace with a strict layered dependency graph. Middle-layer crates must NOT depend on each other — only the top-layer entry points combine them.

```
nexflux-cli / nexflux-python        ← only layer that combines multiple crates
    ├── nexflux-core                ← pure data structures + static reference data
    ├── nexflux-io        → core    ← format readers/writers
    ├── nexflux-structure → core    ← supercell, vacuum, merge, box estimation
    ├── nexflux-analysis  → core    ← md/, dft/ (future), ml/ (future)
    └── nexflux-workflow  → core    ← QC software input file builders
```

### Crate responsibilities

| Crate | Role |
|---|---|
| `nexflux-core` | `Atom`, `Frame`, `Trajectory`, `Cell`; static element/compound data; error types; unit conversion |
| `nexflux-io` | Format readers (`read_xyz`, `read_pdb`, ...) and writers; returns/accepts `Trajectory` |
| `nexflux-structure` | Supercell, vacuum layer, merge, initial box estimation from compound data |
| `nexflux-analysis` | Sub-modules: `md/` (g(r), S(q), MSD, angle, VACF, rotcorr, VanHove, cube density, cluster SDF), `dft/` (future), `ml/` (future) |
| `nexflux-workflow` | QC input file builders: `GaussianJobBuilder`, `GromacsTopologyBuilder`, etc. |
| `nexflux-cli` | CLI + REPL + batch mode (shared interpreter) |
| `nexflux-python` | PyO3 wrappers only; pure Rust libs have zero Python awareness |

---

## Core Data Model

**Primary use case is periodic systems** (crystals, surfaces). Non-periodic (molecular) systems are secondary.

**Always use `Trajectory` as the top-level type**, even for single-frame files.
- Single-frame file → `Trajectory { frames: vec![frame_0] }`
- MD trajectory → `Trajectory { frames: vec![frame_0, frame_1, ...] }`

All module APIs accept/return `Trajectory`. All core types derive `Clone`.

### Atom
```rust
pub struct Atom {
    pub element: String,
    pub position: Vector3<f64>,    // Å, Cartesian
    pub label: Option<String>,     // human-readable tag e.g. "Fe1", "Fe2"
    pub mass: Option<f64>,         // None = look up from elements table
    pub magmom: Option<f64>,       // initial magnetic moment (DFT input)
    pub charge: Option<f64>,       // Bader/DDEC charge (post-processing result)
}
```

Atom index is **implicit** (position in `Vec<Atom>`). No `index` field stored to avoid inconsistency.

### Frame (≈ ASE Atoms)
```rust
pub struct Frame {
    pub atoms: Vec<Atom>,
    pub cell: Option<Cell>,                      // None = non-periodic
    pub pbc: [bool; 3],                          // periodic in x/y/z
    pub charge: i32,                             // total system charge
    pub multiplicity: u32,                       // 2S+1; unpaired electrons = multiplicity-1
    pub bonds: Option<Vec<(usize, usize)>>,      // optional bond list (i,j)
    // post-processing results (written back after calculation):
    pub energy: Option<f64>,
    pub forces: Option<Vec<Vector3<f64>>>,
    pub stress: Option<Matrix3<f64>>,
    pub velocities: Option<Vec<Vector3<f64>>>,
}
```

`pbc` controls periodicity: `[false,false,false]` = molecular, `[true,true,true]` = crystal, `[true,true,false]` = surface/slab.

**`Molecule` struct does not exist** — `Frame` covers all cases.

### nexflux-core file structure
```
nexflux-core/src/
├── lib.rs
├── atom.rs
├── cell.rs
├── frame.rs
├── trajectory.rs
├── data/
│   ├── mod.rs
│   ├── elements.rs  # static: symbol, atomic number, mass, oxidation states, electron config
│   └── compounds.rs # static: name, formula, molecular_mass, density — for box estimation
├── units.rs
└── error.rs
```

### Cell
```rust
pub struct Cell {
    pub matrix: Matrix3<f64>,  // row vectors a, b, c in Å
}
```

| Method | Purpose |
|---|---|
| `from_matrix` / `from_lengths_angles` | constructors |
| `lengths() -> [f64; 3]` | a, b, c |
| `angles() -> [f64; 3]` | α, β, γ in degrees |
| `volume() -> f64` | |
| `fractional_to_cartesian` / `cartesian_to_fractional` | coordinate transforms |
| `wrap_position` | fold Cartesian position back into box |
| `minimum_image` | minimum image convention for PBC distance/vector |

NPT trajectories (varying box per frame) are naturally handled: `cell: Option<Cell>` lives inside `Frame`, so each frame carries its own cell.

### Units
Internal standard follows **DeePMD-kit / VASP convention**:

| Quantity | Internal unit |
|---|---|
| Length | Å |
| Energy | eV |
| Force | eV/Å |
| Stress | eV/Å³ |
| Time | fs |
| Mass | amu |
| Charge | e |
| Temperature | K |

Self-implemented (no `uom` or external unit crates). Enum-based conversion:
```rust
pub enum LengthUnit   { Angstrom, Bohr, Nanometer }
pub enum EnergyUnit   { EV, Hartree, KcalPerMol, KJPerMol, Wavenumber }
pub enum PressureUnit { EVPerAng3, GPa, Kbar }
pub enum TimeUnit     { Femtosecond, Picosecond }

pub fn convert_length(value: f64, from: LengthUnit, to: LengthUnit) -> f64
pub fn convert_energy(value: f64, from: EnergyUnit, to: EnergyUnit) -> f64
pub fn convert_pressure(value: f64, from: PressureUnit, to: PressureUnit) -> f64
```

### Error handling
- Library crates: `nexflux_core::error::ChemError` / `nexflux_core::Result<T>`
- CLI: `anyhow::Result`

---

## Static Reference Data (nexflux-core/data/)

### elements.rs
Per-element: symbol, atomic number, atomic mass, common oxidation states, electron configuration, electronegativity.

### compounds.rs
Used by `nexflux-structure` to estimate initial MD simulation box size:
```rust
pub struct CompoundData {
    pub name: &'static str,
    pub formula: &'static str,
    pub molecular_mass: f64,       // g/mol
    pub density: Option<f64>,      // g/cm³ at standard conditions; None for gases
    pub cas: Option<&'static str>,
}
```
Box estimation logic (V = Σ n_i·M_i / ρ_mix·Nₐ) lives in `nexflux-structure`, not here.

---

## Execution Modes

### 1. One-shot CLI
```bash
nexflux convert -i a.xyz -o b.pdb
```

### 2. Interactive REPL (rustyline)
```bash
nexflux
nexflux> read water.xyz
nexflux> supercell 2 2 1
nexflux> write POSCAR
```

### 3. Batch / script mode
Same interpreter as REPL, non-interactive. For shell script integration.
```bash
nexflux -f workflow.mf
echo -e "read water.xyz\nsupercell 2 2 1\nwrite POSCAR" | nexflux
```

### 4. Python library
`import nexflux` via PyO3 in `nexflux-python`.

### nexflux-cli internal structure
```
nexflux-cli/src/
├── main.rs          # mode detection
├── interpreter.rs   # shared command parser/executor (REPL + batch)
├── repl.rs          # rustyline interactive input
├── batch.rs         # file/stdin input
└── commands/
    ├── io.rs
    ├── structure.rs
    └── analysis.rs
```

---

## Python Bindings (nexflux-python)

All PyO3 glue lives here. Library crates have zero Python awareness.

```
nexflux-python/src/
├── lib.rs        # #[pymodule] entry
├── types.rs      # PyTrajectory (#[pyclass] wrapping inner: Trajectory)
├── io.rs
├── structure.rs
└── analysis.rs
```

Return types: `Vec<f64>`, `HashMap<String, Vec<f64>>` — PyO3 converts automatically to Python list/dict. No numpy or polars Rust crates.

```toml
# nexflux-python/Cargo.toml — minimal deps
[dependencies]
pyo3 = { version = "0.21", features = ["extension-module"] }
nexflux-core      = { path = "../nexflux-core" }
nexflux-io        = { path = "../nexflux-io" }
nexflux-structure = { path = "../nexflux-structure" }
nexflux-analysis  = { path = "../nexflux-analysis" }
```

---

## CLI Binaries Reference

| Binary | Mode / flag | Purpose |
|---|---|---|
| `nex-convert` | `-i <in> -o <out>` | Format conversion |
| `nex-info` | `-i <file>` | Print structure info |
| `nex-job` | `-i <file> -s gaussian [-m B3LYP] [-b 6-31G*]` | Generate QC input file |
| `nex-traj` | `-m gr\|sq\|msd\|angle` | Structural analysis |
| `nex-corr` | `-m vacf\|rotcorr\|vanhove` | Correlation functions |
| `nex-cube` | `-m density\|velocity\|force\|sdf` | Spatial distribution maps |
| `nex-network` | `--P-O=2.3 [--format csv\|xlsx]` | Glass network analysis |

Common flags shared by all trajectory binaries: `--last-n N`, `--ncore N`, `--metal-units`.
Run any binary without `-i` to print mode-specific help.

### `nex-cube -m sdf` — Cluster SDF

Identifies Qn-type clusters (connected components of network-former atoms sharing bridging ligands), aligns each to a reference via Kabsch + permutation enumeration, and outputs per-atom-type 3D probability density as Gaussian cube files.

```
nexflux-analysis/src/md/cube_sdf.rs
```

Key types:
```rust
pub struct ClusterSdfParams {
    pub former: String,              // network former element, e.g. "P"
    pub ligand: String,              // bridging ligand element, e.g. "O"
    pub target_qn: u8,               // 0/1/2/3 — determined by max individual Qn in component
    pub former_ligand_cutoff: f64,
    pub modifier: Option<String>,    // modifier cation, e.g. Some("Zn")
    pub modifier_cutoff: f64,
    pub grid_res: f64,               // Å/voxel
    pub sigma: f64,                  // Gaussian broadening in voxels
    pub padding: f64,                // grid boundary margin in Å
    pub rmsd_warn_threshold: f64,
}
```

Atom-type labels: `P0/P1/P2/P3` (individual Qn), `Of/On/Ob` (O connectivity), modifier element symbol.
Output: `<stem>_<label>.cube` per atom type (multi-family: `<stem>_fam<N>_<label>.cube`).

### `nexflux-analysis/src/md/` file structure

| File | Analysis |
|---|---|
| `gr.rs` | Radial distribution function g(r) + coordination number CN(r) |
| `sq.rs` | Structure factor S(q) via Fourier transform of g(r) |
| `msd.rs` | Mean square displacement (time-shift average, NPT-safe) |
| `angle.rs` | Bond angle distribution P(θ) |
| `vanhove.rs` | Van Hove self-correlation Gs(r, τ) |
| `vacf.rs` | Velocity autocorrelation function |
| `rotcorr.rs` | Rotational correlation C₂(t) for molecular bond vectors |
| `cube_density.rs` | 3D spatial density / velocity / force maps |
| `cube_sdf.rs` | Cluster SDF — Kabsch alignment + per-type density accumulation |

---

## Extending the Project

### Add a file format
1. `nexflux-io/src/readers/<fmt>.rs` — return `Result<Trajectory>`
2. `nexflux-io/src/writers/<fmt>.rs`
3. Export from `readers/mod.rs`, `writers/mod.rs`
4. Add format detection in `nexflux-cli/src/commands/io.rs`
5. Add wrapper in `nexflux-python/src/io.rs`

### Add a structure operation
1. Implement in `nexflux-structure/src/` — takes/returns `Trajectory`
2. `nexflux-cli/src/commands/structure.rs`
3. `nexflux-python/src/structure.rs`

### Add an analysis method
1. Implement in `nexflux-analysis/src/<domain>/`
2. `nexflux-cli/src/commands/analysis.rs`
3. `nexflux-python/src/analysis.rs`

### Add a QC software target
1. Builder in `nexflux-workflow/src/job_builder.rs`
2. Templates in `nexflux-workflow/src/templates.rs`
3. CLI branch in `nexflux-cli/src/commands/`

---

## Key Dependencies

| Crate | Purpose |
|---|---|
| `nalgebra` | 3D vectors/matrices; coordinates are `nalgebra::Vector3<f64>`, cell is `Matrix3<f64>` |
| `ndarray` | Multi-dimensional arrays for bulk trajectory data |
| `rayon` | Data parallelism |
| `thiserror` | Derive macros for `ChemError` |
| `anyhow` | Error propagation in CLI |
| `clap` | CLI argument parsing |
| `rustyline` | Readline-style input for REPL (to be added to workspace deps) |
| `pyo3` | Python bindings (nexflux-python only) |
