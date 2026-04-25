use clap::ValueEnum;
use molflow_analysis::CubeMode;

#[derive(ValueEnum, Clone, Debug)]
pub enum CubeCliMode {
    Density,
    Velocity,
    Force,
}

impl From<CubeCliMode> for CubeMode {
    fn from(m: CubeCliMode) -> Self {
        match m {
            CubeCliMode::Density  => CubeMode::Density,
            CubeCliMode::Velocity => CubeMode::Velocity,
            CubeCliMode::Force    => CubeMode::Force,
        }
    }
}
