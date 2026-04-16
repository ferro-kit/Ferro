//! XYZ 格式文件读取器
//! 
//! XYZ 是最简单的分子结构文件格式，包含原子数、注释行和原子坐标

use molflow_core::{Molecule, Atom};
use nalgebra::Vector3;
use std::fs::File;
use std::io::{BufRead, BufReader};
use anyhow::{Result, Context};

/// 读取 XYZ 格式文件
/// 
/// # 文件格式
/// ```text
/// 3
/// Water molecule
/// O  0.000  0.000  0.000
/// H  0.758  0.587  0.000
/// H -0.758  0.587  0.000
/// ```
/// 
/// # 示例
/// ```no_run
/// use molflow_io::read_xyz;
/// 
/// let molecule = read_xyz("water.xyz").unwrap();
/// println!("Loaded {} atoms", molecule.atom_count());
/// ```
pub fn read_xyz(path: &str) -> Result<Molecule> {
    let file = File::open(path)
        .context(format!("Failed to open file: {}", path))?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    
    // 第一行: 原子数
    let n_atoms: usize = lines.next()
        .context("Empty file")?
        .context("Failed to read first line")?
        .trim()
        .parse()
        .context("Invalid atom count")?;
    
    // 第二行: 注释
    let comment = lines.next()
        .context("Missing comment line")?
        .context("Failed to read comment line")?;
    
    let mut molecule = Molecule::new();
    molecule.set_name(comment.trim());
    
    // 读取原子坐标
    for (i, line) in lines.enumerate() {
        if i >= n_atoms {
            break;
        }
        
        let line = line.context(format!("Failed to read line {}", i + 3))?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        
        if parts.len() < 4 {
            anyhow::bail!("Invalid line format at line {}: expected at least 4 fields", i + 3);
        }
        
        let element = parts[0].to_string();
        let x: f64 = parts[1].parse()
            .context(format!("Invalid x coordinate at line {}", i + 3))?;
        let y: f64 = parts[2].parse()
            .context(format!("Invalid y coordinate at line {}", i + 3))?;
        let z: f64 = parts[3].parse()
            .context(format!("Invalid z coordinate at line {}", i + 3))?;
        
        let position = Vector3::new(x, y, z);
        molecule.add_atom(Atom::new(element, position));
    }
    
    if molecule.atom_count() != n_atoms {
        anyhow::bail!(
            "Atom count mismatch: expected {}, found {}",
            n_atoms,
            molecule.atom_count()
        );
    }
    
    Ok(molecule)
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_xyz_parsing() {
        // 这里可以添加测试，使用临时文件
        // 实际测试需要创建测试文件
    }
}
