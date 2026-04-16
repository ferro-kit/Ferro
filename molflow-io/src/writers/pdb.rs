//! PDB 格式文件写入器

use molflow_core::Molecule;
use std::fs::File;
use std::io::{Write, BufWriter};
use anyhow::{Result, Context};

/// 将分子写入 PDB 格式文件
/// 
/// # 示例
/// ```no_run
/// use molflow_core::Molecule;
/// use molflow_io::write_pdb;
/// 
/// let molecule = Molecule::new();
/// write_pdb(&molecule, "output.pdb").unwrap();
/// ```
pub fn write_pdb(molecule: &Molecule, path: &str) -> Result<()> {
    let file = File::create(path)
        .context(format!("Failed to create file: {}", path))?;
    let mut writer = BufWriter::new(file);
    
    // HEADER 记录
    if let Some(name) = &molecule.name {
        writeln!(writer, "HEADER    {}", name)?;
    }
    
    // ATOM 记录
    for (i, atom) in molecule.atoms.iter().enumerate() {
        // PDB 格式: ATOM   serial atom resName chainID resSeq    x       y       z     occupancy tempFactor element
        writeln!(
            writer,
            "{:<6}{:>5} {:^4} {:3} {:1}{:>4}    {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {:>2}",
            "ATOM",
            i + 1,                    // serial
            atom.element,              // atom name
            "UNK",                    // residue name
            "A",                      // chain ID
            1,                        // residue sequence
            atom.position.x,
            atom.position.y,
            atom.position.z,
            1.00,                     // occupancy
            0.00,                     // temperature factor
            atom.element              // element symbol
        )?;
    }
    
    writeln!(writer, "END")?;
    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_pdb_writing() {
        // 实际测试需要使用临时文件
    }
}
