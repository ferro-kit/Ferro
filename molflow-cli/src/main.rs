fn main() {
    println!(
        r#"molflow — Computational Chemistry Toolkit  v{}

Available commands:
  mol-convert   Format conversion  (xyz, pdb, cif, poscar, extxyz, lammps, qe, ...)
  mol-info      Display structure / trajectory information
  mol-job       Generate QC input files  (gaussian, ...)
  mol-traj      Structural analysis      -m gr | sq | msd | angle
  mol-corr      Correlation functions    -m vacf | rotcorr | vanhove
  mol-cube      Spatial density maps     -m density | velocity | force
  mol-network   Glass network analysis   CN, FO/NBO/BO/OBO, Qn  (--Former-Ligand=cutoff)

Usage:
  <command> --help           Show all flags
  <command> -m <mode>        Show mode-specific help (mol-traj / mol-corr / mol-cube)
  <command> -m <mode> -i ...  Run analysis"#,
        env!("CARGO_PKG_VERSION")
    );
}
