//! 分子数据结构

use crate::Atom;
use nalgebra::Vector3;
use serde::{Serialize, Deserialize};

/// 分子结构
/// 
/// 包含原子列表、键连接信息等
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Molecule {
    /// 原子列表
    pub atoms: Vec<Atom>,
    /// 键连接 (原子索引对)
    pub bonds: Vec<(usize, usize)>,
    /// 分子名称
    pub name: Option<String>,
}

impl Molecule {
    /// 创建空分子
    pub fn new() -> Self {
        Self {
            atoms: Vec::new(),
            bonds: Vec::new(),
            name: None,
        }
    }
    
    /// 添加原子
    pub fn add_atom(&mut self, mut atom: Atom) {
        atom.index = self.atoms.len();
        self.atoms.push(atom);
    }
    
    /// 添加化学键
    pub fn add_bond(&mut self, i: usize, j: usize) {
        if i < self.atoms.len() && j < self.atoms.len() {
            self.bonds.push((i, j));
        }
    }
    
    /// 设置分子名称
    pub fn set_name(&mut self, name: impl Into<String>) {
        self.name = Some(name.into());
    }
    
    /// 获取原子数量
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }
    
    /// 获取键数量
    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }
    
    /// 计算质心
    pub fn center_of_mass(&self) -> Vector3<f64> {
        if self.atoms.is_empty() {
            return Vector3::zeros();
        }
        
        let total_mass: f64 = self.atoms.iter().map(|a| a.mass).sum();
        let weighted_sum = self.atoms.iter()
            .map(|a| a.position * a.mass)
            .fold(Vector3::zeros(), |acc, v| acc + v);
        
        weighted_sum / total_mass
    }
    
    /// 计算几何中心
    pub fn geometric_center(&self) -> Vector3<f64> {
        if self.atoms.is_empty() {
            return Vector3::zeros();
        }
        
        let sum = self.atoms.iter()
            .map(|a| a.position)
            .fold(Vector3::zeros(), |acc, v| acc + v);
        
        sum / (self.atoms.len() as f64)
    }
    
    /// 平移分子
    pub fn translate(&mut self, displacement: Vector3<f64>) {
        for atom in &mut self.atoms {
            atom.position += displacement;
        }
    }
    
    /// 将分子质心移到原点
    pub fn center(&mut self) {
        let com = self.center_of_mass();
        self.translate(-com);
    }
    
    /// 计算总质量
    pub fn total_mass(&self) -> f64 {
        self.atoms.iter().map(|a| a.mass).sum()
    }
    
    /// 获取指定元素的原子数量
    pub fn count_element(&self, element: &str) -> usize {
        self.atoms.iter().filter(|a| a.element == element).count()
    }
}

impl Default for Molecule {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_molecule_creation() {
        let mol = Molecule::new();
        assert_eq!(mol.atom_count(), 0);
        assert_eq!(mol.bond_count(), 0);
    }
    
    #[test]
    fn test_add_atoms() {
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new("C", Vector3::new(0.0, 0.0, 0.0)));
        mol.add_atom(Atom::new("H", Vector3::new(1.0, 0.0, 0.0)));
        assert_eq!(mol.atom_count(), 2);
    }
    
    #[test]
    fn test_center_of_mass() {
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new("C", Vector3::new(0.0, 0.0, 0.0)));
        mol.add_atom(Atom::new("C", Vector3::new(2.0, 0.0, 0.0)));
        let com = mol.center_of_mass();
        assert!((com.x - 1.0).abs() < 1e-10);
    }
}
