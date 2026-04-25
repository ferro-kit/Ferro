//! nexflux-core — 核心数据结构与工具
//!
//! 层次结构：
//! ```text
//! Trajectory (Vec<Frame>)
//!   └── Frame  (atoms + cell + pbc + 计算结果)
//!         └── Atom (element + position + 可选属性)
//!               └── data::elements (静态元素表，用于质量查找)
//! ```

pub mod atom;
pub mod cell;
pub mod frame;
pub mod trajectory;
pub mod cube_data;
pub mod data;
pub mod units;
pub mod error;

// 顶层重导出，方便下游 crate 使用
pub use atom::Atom;
pub use cell::Cell;
pub use frame::Frame;
pub use trajectory::{Trajectory, TrajectoryMetadata};
pub use cube_data::CubeData;
pub use error::{ChemError, Result};
