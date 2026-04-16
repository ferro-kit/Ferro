//! 几何分析工具

use molflow_core::Molecule;
use nalgebra::Vector3;

/// 计算两个原子之间的距离
pub fn distance(mol: &Molecule, i: usize, j: usize) -> Option<f64> {
    let atom1 = mol.atoms.get(i)?;
    let atom2 = mol.atoms.get(j)?;
    Some((atom1.position - atom2.position).norm())
}

/// 计算三个原子之间的键角 (返回角度，单位：度)
pub fn angle(mol: &Molecule, i: usize, j: usize, k: usize) -> Option<f64> {
    let atom1 = mol.atoms.get(i)?;
    let atom2 = mol.atoms.get(j)?;
    let atom3 = mol.atoms.get(k)?;
    
    let v1 = atom1.position - atom2.position;
    let v2 = atom3.position - atom2.position;
    
    let cos_angle = v1.dot(&v2) / (v1.norm() * v2.norm());
    Some(cos_angle.acos().to_degrees())
}

/// 计算四个原子之间的二面角 (返回角度，单位：度)
pub fn dihedral(mol: &Molecule, i: usize, j: usize, k: usize, l: usize) -> Option<f64> {
    let atom1 = mol.atoms.get(i)?;
    let atom2 = mol.atoms.get(j)?;
    let atom3 = mol.atoms.get(k)?;
    let atom4 = mol.atoms.get(l)?;
    
    let b1 = atom2.position - atom1.position;
    let b2 = atom3.position - atom2.position;
    let b3 = atom4.position - atom3.position;
    
    let n1 = b1.cross(&b2);
    let n2 = b2.cross(&b3);
    
    let m1 = n1.cross(&b2.normalize());
    
    let x = n1.dot(&n2);
    let y = m1.dot(&n2);
    
    Some(y.atan2(x).to_degrees())
}

/// 计算分子的回转半径
pub fn radius_of_gyration(mol: &Molecule) -> f64 {
    if mol.atoms.is_empty() {
        return 0.0;
    }
    
    let com = mol.center_of_mass();
    let total_mass = mol.total_mass();
    
    let sum: f64 = mol.atoms.iter()
        .map(|atom| {
            let r = (atom.position - com).norm();
            atom.mass * r * r
        })
        .sum();
    
    (sum / total_mass).sqrt()
}

/// 计算最小包围盒尺寸
pub fn bounding_box(mol: &Molecule) -> Option<Vector3<f64>> {
    if mol.atoms.is_empty() {
        return None;
    }
    
    let first = &mol.atoms[0].position;
    let mut min = *first;
    let mut max = *first;
    
    for atom in &mol.atoms {
        min.x = min.x.min(atom.position.x);
        min.y = min.y.min(atom.position.y);
        min.z = min.z.min(atom.position.z);
        
        max.x = max.x.max(atom.position.x);
        max.y = max.y.max(atom.position.y);
        max.z = max.z.max(atom.position.z);
    }
    
    Some(max - min)
}

#[cfg(test)]
mod tests {
    use super::*;
    use molflow_core::Atom;
    
    #[test]
    fn test_distance() {
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new("C", Vector3::new(0.0, 0.0, 0.0)));
        mol.add_atom(Atom::new("C", Vector3::new(1.0, 0.0, 0.0)));
        
        let d = distance(&mol, 0, 1).unwrap();
        assert!((d - 1.0).abs() < 1e-10);
    }
}
