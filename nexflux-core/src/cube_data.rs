//! CubeData — Gaussian cube 文件的完整数据载体
//!
//! Cube 文件包含两部分：原子结构（单帧）和三维体积数据（电子密度、波函数等）。
//! 两者不可分离，故用独立结构体承载，不污染通用的 [`Frame`]。

use ndarray::Array3;
use nalgebra::{Matrix3, Vector3};
use crate::Frame;

/// Container for a Gaussian cube file: atomic structure plus volumetric grid data.
///
/// All length quantities are in **Ångström** (internal standard units).
///
/// The volumetric data is stored in C-order (row-major) with axis order X → Y → Z,
/// matching the cube format convention (outer loop X, middle loop Y, inner loop Z).
#[derive(Debug, Clone)]
pub struct CubeData {
    /// Atomic structure: atoms, cell, pbc (single frame)
    pub frame: Frame,
    /// Volumetric data array, shape `(nx, ny, nz)` \[arbitrary units, e.g. e/Å³\]
    pub data: Array3<f64>,
    /// Grid origin \[Å\]
    pub origin: Vector3<f64>,
    /// Voxel vectors \[Å\]: row `i` is the step vector along grid axis `i`
    pub spacing: Matrix3<f64>,
}

impl CubeData {
    /// Grid dimensions `(nx, ny, nz)`.
    pub fn shape(&self) -> (usize, usize, usize) {
        let s = self.data.shape();
        (s[0], s[1], s[2])
    }
}
