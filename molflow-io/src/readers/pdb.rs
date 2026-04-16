//! PDB (Protein Data Bank) 格式文件读取器
//! 
//! PDB 是生物大分子结构的标准格式

use molflow_core::{Molecule, Atom};
use nalgebra::Vector3;
use std::fs::File;
use std::io::{BufRead, BufReader};
use anyhow::{Result, Context};

/// 读取 PDB 格式文件
/// 
/// # 文件格式示例
/// ```text
/// ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O
/// ATOM      2  H1  HOH A   1       0.758   0.587   0.000  1.00  0.00           H
/// ATOM      3  H2  HOH A   1      -0.758   0.587   0.000  1.00  0.00           H
/// ```
/// 
/// # 示例
/// ```no_run
/// use molflow_io::read_pdb;
/// 
/// let molecule = read_pdb("protein.pdb").unwrap();
/// println!("Loaded {} atoms", molecule.atom_count());
/// ```
pub fn read_pdb(path: &str) -> Result<Molecule> {
    let file = File::open(path)
        .context(format!("Failed to open file: {}", path))?;
    let reader = BufReader::new(file);
    
    let mut molecule = Molecule::new();
    
    for line in reader.lines() {
        let line = line?;
        
        // 只处理 ATOM 和 HETATM 记录
        if !line.starts_with("ATOM") && !line.starts_with("HETATM") {
            continue;
        }
        
        // PDB 是固定宽度格式
        if line.len() < 54 {
            continue;
        }
        
        // 提取元素符号 (列 77-78，如果存在)
        let element = if line.len() >= 78 {
            line[76..78].trim().to_string()
        } else {
            // 尝试从原子名称推断元素
            line[12..16].trim().chars().next()
                .map(|c| c.to_string())
                .unwrap_or_else(|| "X".to_string())
        };
        
        // 提取坐标 (列 31-54)
        let x: f64 = line[30..38].trim().parse()
            .context("Invalid x coordinate")?;
        let y: f64 = line[38..46].trim().parse()
            .context("Invalid y coordinate")?;
        let z: f64 = line[46..54].trim().parse()
            .context("Invalid z coordinate")?;
        
        let position = Vector3::new(x, y, z);
        molecule.add_atom(Atom::new(element, position));
    }
    
    Ok(molecule)
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_pdb_parsing() {
        // 实际测试需要创建测试文件
    }
}
