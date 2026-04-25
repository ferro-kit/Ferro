use clap::ValueEnum;
use nexflux_analysis::CubeMode;

#[derive(ValueEnum, Clone, Debug)]
pub enum CubeCliMode {
    Density,
    Velocity,
    Force,
    Sdf,
}

impl From<CubeCliMode> for CubeMode {
    fn from(m: CubeCliMode) -> Self {
        match m {
            CubeCliMode::Density  => CubeMode::Density,
            CubeCliMode::Velocity => CubeMode::Velocity,
            CubeCliMode::Force    => CubeMode::Force,
            CubeCliMode::Sdf      => CubeMode::Density, // unreachable; sdf takes a separate path
        }
    }
}
