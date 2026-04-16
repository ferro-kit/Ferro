//! ChemTools Core - 核心数据结构和类型定义
//! 
//! 这个模块提供了计算化学中的基本数据结构，如原子、分子、轨迹等。

pub mod atom;
pub mod molecule;
pub mod trajectory;
pub mod units;
pub mod error;

// 重新导出常用类型，方便使用
pub use atom::Atom;
pub use molecule::Molecule;
pub use trajectory::Trajectory;
pub use error::{ChemError, Result};

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    #[test]
    fn test_basic_molecule() {
        let mut mol = Molecule::new();
        let atom1 = Atom::new("C", Vector3::new(0.0, 0.0, 0.0));
        let atom2 = Atom::new("H", Vector3::new(1.0, 0.0, 0.0));
        
        mol.add_atom(atom1);
        mol.add_atom(atom2);
        
        assert_eq!(mol.atom_count(), 2);
    }
}
