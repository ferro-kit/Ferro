use anyhow::{anyhow, Result};
use clap::Parser;
use molflow::{
    args::cube::CubeCliMode,
    help::print_cube_help,
    io_dispatch::read_trajectory,
};
use molflow_analysis::{calc_cube_density, CubeDensityParams};
use molflow_io::{write_cube, LammpsUnits};
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[command(
    name = "mol-cube",
    about = "Spatial distribution maps: density | velocity | force  (run without -i for mode help)"
)]
struct Cli {
    /// Analysis mode
    #[arg(short = 'm', long, value_enum)]
    mode: CubeCliMode,

    /// Input trajectory file (omit to show mode-specific help)
    #[arg(short, long)]
    input: Option<PathBuf>,

    /// Output cube file
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Use only the last N frames
    #[arg(long)]
    last_n: Option<usize>,

    /// Parallel threads (default: all cores)
    #[arg(long)]
    ncore: Option<usize>,

    /// Grid points along a axis
    #[arg(long, default_value = "50")]
    nx: usize,

    /// Grid points along b axis
    #[arg(long, default_value = "50")]
    ny: usize,

    /// Grid points along c axis
    #[arg(long, default_value = "50")]
    nz: usize,

    /// Only include these elements, e.g. Fe,O
    #[arg(long, value_delimiter = ',')]
    elements: Option<Vec<String>>,

    /// Use LAMMPS metal units for dump files (velocities Å/ps, forces eV/Å)
    #[arg(long)]
    metal_units: bool,
}

fn main() -> Result<()> {
    let args = Cli::parse();

    let input = match &args.input {
        Some(p) => p.clone(),
        None => {
            print_cube_help(&args.mode);
            return Ok(());
        }
    };

    if let Some(n) = args.ncore {
        rayon::ThreadPoolBuilder::new().num_threads(n).build_global().ok();
    }

    let units = if args.metal_units { LammpsUnits::Metal } else { LammpsUnits::Real };
    let mut traj = read_trajectory(&input, units)?;
    if let Some(n) = args.last_n {
        traj = traj.tail(n);
    }

    let default_name = format!("{}.cube", mode_name(&args.mode));
    let out = args.output.as_deref().unwrap_or(Path::new(&default_name));

    let params = CubeDensityParams {
        nx: args.nx,
        ny: args.ny,
        nz: args.nz,
        elements: args.elements.clone(),
        mode: args.mode.clone().into(),
    };

    let result = calc_cube_density(&traj, &params)
        .ok_or_else(|| anyhow!("Cube calc failed (missing cell, velocities, or forces?)"))?;

    write_cube(out.to_str().unwrap_or(&default_name), &result.cube)?;
    println!(
        "Cube ({}) -> {}  [{} frames, {} atoms]",
        mode_name(&args.mode),
        out.display(),
        result.n_frames,
        result.n_atoms,
    );
    Ok(())
}

fn mode_name(m: &CubeCliMode) -> &'static str {
    match m {
        CubeCliMode::Density  => "density",
        CubeCliMode::Velocity => "velocity",
        CubeCliMode::Force    => "force",
    }
}
