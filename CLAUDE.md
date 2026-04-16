# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Test Commands

```bash
cargo build                                                    # build entire workspace
cargo build --release
cargo build --package molflow-core                             # single crate
cargo test                                                     # all tests
cargo test --package molflow-io                                # single crate
cargo test --package molflow-core test_basic_molecule          # single test
cargo fmt
cargo clippy
cargo check

# CLI
cargo run --bin molflow -- info -i examples/water.xyz
cargo run --bin molflow -- convert -i input.xyz -o output.pdb
cargo run --bin molflow -- job -i input.xyz -s gaussian -m B3LYP -o job.gjf

# Python bindings (requires maturin; molflow-python not yet in workspace)
cd molflow-python && maturin develop
```

Test fixtures: `tests/` (water.xyz, water.pdb, water.cif)

---

## Architecture

Cargo workspace with a strict layered dependency graph. Middle-layer crates must NOT depend on each other — only the top-layer entry points combine them.

```
molflow-cli / molflow-python        ← only layer that combines multiple crates
    ├── molflow-core                ← pure data structures + static reference data
    ├── molflow-io        → core    ← format readers/writers
    ├── molflow-structure → core    ← supercell, vacuum, merge, box estimation
    ├── molflow-analysis  → core    ← md/, dft/ (future), ml/ (future)
    └── molflow-workflow  → core    ← QC software input file builders
```

### Crate responsibilities

| Crate | Role |
|---|---|
| `molflow-core` | `Atom`, `Frame`, `Trajectory`, `Cell`; static element/compound data; error types; unit conversion |
| `molflow-io` | Format readers (`read_xyz`, `read_pdb`, ...) and writers; returns/accepts `Trajectory` |
| `molflow-structure` | Supercell, vacuum layer, merge, initial box estimation from compound data |
| `molflow-analysis` | Sub-modules: `md/` (RMSD, MSD, Rg, RDF), `dft/` (future), `ml/` (future) |
| `molflow-workflow` | QC input file builders: `GaussianJobBuilder`, `GromacsTopologyBuilder`, etc. |
| `molflow-cli` | CLI + REPL + batch mode (shared interpreter) |
| `molflow-python` | PyO3 wrappers only; pure Rust libs have zero Python awareness |

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

### molflow-core file structure
```
molflow-core/src/
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
- Library crates: `molflow_core::error::ChemError` / `molflow_core::Result<T>`
- CLI: `anyhow::Result`

---

## Static Reference Data (molflow-core/data/)

### elements.rs
Per-element: symbol, atomic number, atomic mass, common oxidation states, electron configuration, electronegativity.

### compounds.rs
Used by `molflow-structure` to estimate initial MD simulation box size:
```rust
pub struct CompoundData {
    pub name: &'static str,
    pub formula: &'static str,
    pub molecular_mass: f64,       // g/mol
    pub density: Option<f64>,      // g/cm³ at standard conditions; None for gases
    pub cas: Option<&'static str>,
}
```
Box estimation logic (V = Σ n_i·M_i / ρ_mix·Nₐ) lives in `molflow-structure`, not here.

---

## Execution Modes

### 1. One-shot CLI
```bash
molflow convert -i a.xyz -o b.pdb
```

### 2. Interactive REPL (rustyline)
```bash
molflow
molflow> read water.xyz
molflow> supercell 2 2 1
molflow> write POSCAR
```

### 3. Batch / script mode
Same interpreter as REPL, non-interactive. For shell script integration.
```bash
molflow -f workflow.mf
echo -e "read water.xyz\nsupercell 2 2 1\nwrite POSCAR" | molflow
```

### 4. Python library
`import molflow` via PyO3 in `molflow-python`.

### molflow-cli internal structure
```
molflow-cli/src/
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

## Python Bindings (molflow-python)

All PyO3 glue lives here. Library crates have zero Python awareness.

```
molflow-python/src/
├── lib.rs        # #[pymodule] entry
├── types.rs      # PyTrajectory (#[pyclass] wrapping inner: Trajectory)
├── io.rs
├── structure.rs
└── analysis.rs
```

Return types: `Vec<f64>`, `HashMap<String, Vec<f64>>` — PyO3 converts automatically to Python list/dict. No numpy or polars Rust crates.

```toml
# molflow-python/Cargo.toml — minimal deps
[dependencies]
pyo3 = { version = "0.21", features = ["extension-module"] }
molflow-core      = { path = "../molflow-core" }
molflow-io        = { path = "../molflow-io" }
molflow-structure = { path = "../molflow-structure" }
molflow-analysis  = { path = "../molflow-analysis" }
```

---

## Extending the Project

### Add a file format
1. `molflow-io/src/readers/<fmt>.rs` — return `Result<Trajectory>`
2. `molflow-io/src/writers/<fmt>.rs`
3. Export from `readers/mod.rs`, `writers/mod.rs`
4. Add format detection in `molflow-cli/src/commands/io.rs`
5. Add wrapper in `molflow-python/src/io.rs`

### Add a structure operation
1. Implement in `molflow-structure/src/` — takes/returns `Trajectory`
2. `molflow-cli/src/commands/structure.rs`
3. `molflow-python/src/structure.rs`

### Add an analysis method
1. Implement in `molflow-analysis/src/<domain>/`
2. `molflow-cli/src/commands/analysis.rs`
3. `molflow-python/src/analysis.rs`

### Add a QC software target
1. Builder in `molflow-workflow/src/job_builder.rs`
2. Templates in `molflow-workflow/src/templates.rs`
3. CLI branch in `molflow-cli/src/commands/`

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
| `pyo3` | Python bindings (molflow-python only) |
