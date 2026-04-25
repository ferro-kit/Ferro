//! MD 轨迹后处理分析子模块
//!
//! 当前实现：
//!   - [`gr`]       — 径向分布函数 g(r) 及配位数 CN(r)
//!   - [`sq`]       — 结构因子 S(q)（由 g(r) Fourier 变换得到）
//!   - [`msd`]      — 均方位移 MSD（time-shift 平均，支持 NPT）
//!   - [`angle`]    — 键角分布
//!   - [`vanhove`]  — van Hove 自相关函数 Gs(r, τ)
//!
pub mod gr;
pub mod sq;
pub mod msd;
pub mod angle;
pub mod vanhove;
pub mod vacf;
pub mod rotcorr;
pub mod cube_density;
pub mod cube_radius;
pub mod cube_sdf;
pub mod scattering_data;

pub use gr::{GrParams, GrResult, PairStats, calc_gr, write_gr, write_cn};
pub use sq::{SqParams, SqResult, SqWeighting, calc_sq_from_gr, write_sq};
pub use msd::{MsdParams, MsdResult, calc_msd, write_msd};
pub use angle::{AngleParams, AngleResult, AngleStats, calc_angle, write_angle};
pub use vanhove::{VanHoveParams, VanHoveResult, calc_vanhove, write_vanhove};
pub use vacf::{VacfParams, VacfResult, calc_vacf, write_vacf};
pub use rotcorr::{RotCorrParams, RotCorrResult, calc_rotcorr, write_rotcorr};
pub use cube_density::{CubeMode, CubeDensityParams, CubeDensityResult, calc_cube_density};
pub use cube_radius::{CubeRadiusParams, CubeRadiusResult, calc_cube_radius};
pub use cube_sdf::{ClusterSdfParams, ClusterFamily, RmsdStats, ClusterSdfResult, calc_cluster_sdf};
