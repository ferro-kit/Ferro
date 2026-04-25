fn main() {
    println!(
        r#"nexflux — Computational Chemistry Toolkit  v{}

Available commands:
  nex-convert   Format conversion  (xyz, pdb, cif, poscar, extxyz, lammps, qe, ...)
  nex-info      Display structure / trajectory information
  nex-job       Generate QC input files  (gaussian, ...)
  nex-traj      Structural analysis      -m gr | sq | msd | angle
  nex-corr      Correlation functions    -m vacf | rotcorr | vanhove
  nex-cube      Spatial density maps     -m density | velocity | force
  nex-network   Glass network analysis   CN, FO/NBO/BO/OBO, Qn  (--Former-Ligand=cutoff)

Usage:
  <command> --help           Show all flags
  <command> -m <mode>        Show mode-specific help (nex-traj / nex-corr / nex-cube)
  <command> -m <mode> -i ...  Run analysis"#,
        env!("CARGO_PKG_VERSION")
    );
}
