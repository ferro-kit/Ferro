//! 计算任务构建器
//! 
//! 用于创建各种量子化学和分子动力学软件的输入文件

use molflow_core::Molecule;
use anyhow::Result;

/// Gaussian 输入文件构建器
pub struct GaussianJobBuilder {
    pub molecule: Molecule,
    pub method: String,
    pub basis_set: String,
    pub charge: i32,
    pub multiplicity: i32,
}

impl GaussianJobBuilder {
    pub fn new(molecule: Molecule) -> Self {
        Self {
            molecule,
            method: "B3LYP".to_string(),
            basis_set: "6-31G*".to_string(),
            charge: 0,
            multiplicity: 1,
        }
    }
    
    /// 生成 Gaussian 输入文件
    pub fn build(&self) -> Result<String> {
        let mut output = String::new();
        output.push_str(&format!("%chk=job.chk\n"));
        output.push_str(&format!("# {} {}\n\n", self.method, self.basis_set));
        output.push_str("Title\n\n");
        output.push_str(&format!("{} {}\n", self.charge, self.multiplicity));
        
        for atom in &self.molecule.atoms {
            output.push_str(&format!(
                "{:<2}  {:>12.6}  {:>12.6}  {:>12.6}\n",
                atom.element, atom.position.x, atom.position.y, atom.position.z
            ));
        }
        output.push_str("\n");
        
        Ok(output)
    }
}

/// GROMACS 拓扑文件构建器 (示例)
pub struct GromacsTopologyBuilder {
    pub molecule: Molecule,
}

impl GromacsTopologyBuilder {
    pub fn new(molecule: Molecule) -> Self {
        Self { molecule }
    }
    
    /// 生成 GROMACS 拓扑文件 (简化版)
    pub fn build(&self) -> Result<String> {
        // TODO: 实现完整的拓扑生成
        Ok("[ moleculetype ]\n; Name  nrexcl\nMOL     3\n".to_string())
    }
}
