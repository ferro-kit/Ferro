# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Test Commands

```bash
# Build entire workspace
cargo build
cargo build --release

# Build a single crate
cargo build --package molflow-core

# Run all tests
cargo test

# Run tests for a specific crate
cargo test --package molflow-io

# Run a single test by name
cargo test --package molflow-core test_basic_molecule

# Format code
cargo fmt

# Lint / static analysis
cargo clippy

# Type-check without compiling
cargo check

# Run the CLI
cargo run --bin molflow -- info -i examples/water.xyz
cargo run --bin molflow -- convert -i input.xyz -o output.pdb
cargo run --bin molflow -- job -i input.xyz -s gaussian -m B3LYP -o job.gjf

# Python bindings (optional, requires maturin)
cd molflow-python && maturin develop
```

Test fixtures are in `tests/` (water.xyz, water.pdb, water.cif).

## Architecture

This is a Cargo **workspace** of 5 active crates (plus a commented-out `molflow-python`). The dependency graph flows strictly downward:

```
molflow-cli
  ├── molflow-core      (lowest layer — no internal deps)
  ├── molflow-io        → molflow-core
  ├── molflow-analysis  → molflow-core
  └── molflow-workflow  → molflow-core
```

**Do not introduce cross-dependencies between the middle layers** (e.g., `molflow-analysis` must not depend on `molflow-io`).

### Crate responsibilities

| Crate | Role |
|---|---|
| `molflow-core` | `Atom`, `Molecule`, `Trajectory` structs; `ChemError`/`Result` types; unit conversion |
| `molflow-io` | File format readers (`read_xyz`, `read_pdb`) and writers (`write_xyz`, `write_pdb`) |
| `molflow-analysis` | Geometry (bond lengths/angles/dihedrals), trajectory (RMSD, MSD, Rg), properties (dipole) |
| `molflow-workflow` | Input file builders for QC software: `GaussianJobBuilder`, `GromacsTopologyBuilder` |
| `molflow-cli` | `clap`-based CLI with subcommands: `convert`, `analyze`, `job`, `info` |

### Error handling convention

- **Library crates** (`core`, `io`, `analysis`, `workflow`) return `molflow_core::error::ChemError` / `molflow_core::Result<T>`.
- **CLI** (`molflow-cli`) uses `anyhow::Result` for ergonomic error propagation.

### Adding a new file format

1. `molflow-io/src/readers/<format>.rs` — implement reader, return `molflow_core::Result<Molecule>`
2. `molflow-io/src/writers/<format>.rs` — implement writer
3. Export from `molflow-io/src/readers/mod.rs` and `writers/mod.rs`
4. Add format detection in `molflow-cli/src/main.rs`

### Adding a new QC software target

1. Add a builder struct in `molflow-workflow/src/job_builder.rs`
2. Register templates in `molflow-workflow/src/templates.rs`
3. Add a CLI branch in `create_job()` in `molflow-cli/src/main.rs`

## Key Dependencies

| Crate | Purpose |
|---|---|
| `nalgebra` | 3D vectors/matrices (coordinates are `nalgebra::Vector3<f64>`) |
| `ndarray` | Multi-dimensional arrays for bulk data |
| `rayon` | Data parallelism |
| `thiserror` | Derive macros for `ChemError` |
| `anyhow` | Error propagation in CLI |
| `clap` | CLI argument parsing |
| `pyo3` | Python bindings (molflow-python, not in active workspace) |
