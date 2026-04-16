//! 原子数据结构

use nalgebra::Vector3;
use serde::{Serialize, Deserialize};

/// 原子结构
/// 
/// 存储原子的基本信息：元素类型、位置、速度、电荷、质量等
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Atom {
    /// 元素符号 (如 "C", "H", "O")
    pub element: String,
    /// 原子位置坐标 (单位: Å)
    pub position: Vector3<f64>,
    /// 原子速度 (单位: Å/ps) - 可选
    pub velocity: Option<Vector3<f64>>,
    /// 原子电荷 (单位: e) - 可选
    pub charge: Option<f64>,
    /// 原子质量 (单位: amu)
    pub mass: f64,
    /// 原子索引
    pub index: usize,
}

impl Atom {
    /// 创建新原子
    /// 
    /// # 示例
    /// ```
    /// use molflow_core::Atom;
    /// use nalgebra::Vector3;
    /// 
    /// let atom = Atom::new("C", Vector3::new(0.0, 0.0, 0.0));
    /// ```
    pub fn new(element: impl Into<String>, position: Vector3<f64>) -> Self {
        let element = element.into();
        let mass = Self::element_mass(&element);
        Self {
            element,
            position,
            velocity: None,
            charge: None,
            mass,
            index: 0,
        }
    }
    
    /// 根据元素符号获取标准原子质量
    fn element_mass(element: &str) -> f64 {
        match element {
            "H" => 1.008,
            "C" => 12.011,
            "N" => 14.007,
            "O" => 15.999,
            "P" => 30.974,
            "S" => 32.06,
            "F" => 18.998,
            "Cl" => 35.45,
            "Br" => 79.904,
            "I" => 126.90,
            "Na" => 22.990,
            "K" => 39.098,
            "Ca" => 40.078,
            "Mg" => 24.305,
            "Fe" => 55.845,
            "Cu" => 63.546,
            "Zn" => 65.38,
            _ => 1.0, // 默认质量
        }
    }
    
    /// 计算与另一个原子的距离
    pub fn distance_to(&self, other: &Atom) -> f64 {
        (self.position - other.position).norm()
    }
    
    /// 设置原子速度
    pub fn set_velocity(&mut self, velocity: Vector3<f64>) {
        self.velocity = Some(velocity);
    }
    
    /// 设置原子电荷
    pub fn set_charge(&mut self, charge: f64) {
        self.charge = Some(charge);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_atom_creation() {
        let atom = Atom::new("C", Vector3::new(0.0, 0.0, 0.0));
        assert_eq!(atom.element, "C");
        assert_eq!(atom.mass, 12.011);
    }
    
    #[test]
    fn test_distance() {
        let atom1 = Atom::new("C", Vector3::new(0.0, 0.0, 0.0));
        let atom2 = Atom::new("H", Vector3::new(1.0, 0.0, 0.0));
        assert!((atom1.distance_to(&atom2) - 1.0).abs() < 1e-10);
    }
}
