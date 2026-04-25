# nexflux

A modular computational chemistry toolkit written in Rust, designed primarily for **periodic systems** (crystals, surfaces, glasses).

## Workspace Layout

```
nexflux/
├── nexflux-core/       # Core data structures (Atom, Frame, Trajectory, Cell)
├── nexflux-io/         # File readers and writers (10+ formats)
├── nexflux-analysis/   # Post-processing (MD analysis, glass network analysis)
├── nexflux-workflow/   # QC input file generation (Gaussian, etc.)
├── nexflux-cli/        # Command-line binaries
└── nexflux-python/     # Python bindings via PyO3 (not yet enabled)
```

Dependency graph — strictly one-directional; middle-layer crates must not depend on each other:

```
nexflux-cli / nexflux-python
    ├── nexflux-core
    ├── nexflux-io        → core
    ├── nexflux-analysis  → core
    └── nexflux-workflow  → core
```

## Build

```bash
cargo build --release
```

Binaries are placed in `target/release/`:
`nex-convert` · `nex-info` · `nex-job` · `nex-traj` · `nex-corr` · `nex-cube` · `nex-network`

## Supported File Formats

| Format | Read | Write |
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

## CLI Reference

### nex-convert — Format conversion
```bash
nex-convert -i structure.xyz -o structure.pdb
nex-convert -i traj.dump -o traj.extxyz
```

### nex-info — Print structure information
```bash
nex-info -i water.xyz
nex-info -i NaCl.cif
```

### nex-job — Generate QC input files
```bash
nex-job -i water.xyz -s gaussian -m B3LYP -b "6-31G*" -o job.gjf
```

### nex-traj — Structural trajectory analysis
```bash
nex-traj -m gr    -i traj.dump                         # Radial distribution function g(r)
nex-traj -m sq    -i traj.dump --weighting both        # Structure factor S(q)
nex-traj -m msd   -i traj.dump --dt 2.0 --elements Li  # Mean square displacement
nex-traj -m angle -i traj.dump --r-cut-ab 2.0          # Bond-angle distribution
nex-traj -m gr                                         # Run without -i to show mode help
```

### nex-corr — Correlation functions
```bash
nex-corr -m vacf    -i traj.dump --dt 2.0               # Velocity autocorrelation
nex-corr -m rotcorr -i traj.xyz --center O --neighbor H  # Rotational correlation C₂(t)
nex-corr -m vanhove -i traj.dump --tau 100               # Van Hove self-correlation Gs(r,τ)
```

### nex-cube — Spatial density maps
```bash
nex-cube -m density  -i traj.dump                # Number density (atoms/Å³)
nex-cube -m velocity -i traj.dump                # Velocity field
nex-cube -m force    -i traj.dump --elements O   # Force field
nex-cube -m sdf      -i traj.dump --qn 3 --modifier Zn  # Cluster SDF (Kabsch-aligned)
```

Output is Gaussian cube format, readable by VESTA and VMD.

### nex-network — Glass network analysis

Analyzes the connectivity topology between network formers and ligands in periodic systems (glasses, crystals).

```bash
# P–O system (cutoff 2.3 Å)
nex-network -i traj.dump --P-O=2.3

# Multi-element system, xlsx output
nex-network -i traj.dump --P-O=2.3 --Si-O=1.8 --Al-O=2.0 --format xlsx -o result

# Pair flag format: --Former-Ligand=cutoff_in_Angstrom  (Former is capitalized)
```

Output columns:

| Output | Description |
|---|---|
| CN distribution | Coordination number per (former, ligand) pair + mean |
| Ligand classes | FO / NBO(X) / BO(X-Y) / OBO(X-Y-Z) per ligand atom |
| Qn species | Q^n(mX) distribution — n = bridging ligands, m = hetero-element bridges |

CSV: three files (`_cn.csv`, `_ligand.csv`, `_qn.csv`). Excel: three sheets in one `.xlsx`.

## Internal Units

| Quantity | Unit |
|---|---|
| Length | Å |
| Energy | eV |
| Force | eV/Å |
| Stress | eV/Å³ |
| Time | fs |
| Mass | amu |
| Charge | e |

## Core Data Model

```rust
// Trajectory is always the top-level type, even for single-frame files.
pub struct Trajectory {
    pub frames: Vec<Frame>,
    pub metadata: TrajectoryMetadata,
}

pub struct Frame {
    pub atoms: Vec<Atom>,
    pub cell: Option<Cell>,                 // None = non-periodic
    pub pbc: [bool; 3],                     // PBC per axis
    pub energy: Option<f64>,
    pub forces: Option<Vec<Vector3<f64>>>,
    pub velocities: Option<Vec<Vector3<f64>>>,
    // ...
}

pub struct Atom {
    pub element: String,
    pub position: Vector3<f64>,             // Å, Cartesian
    pub label: Option<String>,
    pub mass: Option<f64>,
    pub magmom: Option<f64>,
    pub charge: Option<f64>,
}
```

## Development Commands

```bash
cargo build                           # Build workspace
cargo build --release                 # Release build
cargo test                            # Run all tests
cargo test --package nexflux-io       # Test a single crate
cargo fmt && cargo clippy             # Format and lint
cargo run --bin nex-info -- -i examples/water.xyz
```

## Example Files

`examples/` contains sample inputs: LAMMPS dump, CIF, LAMMPS data, and CP2K `.inp` files.
