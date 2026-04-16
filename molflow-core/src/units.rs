//! 单位转换工具

/// 长度单位转换
pub mod length {
    /// Angstrom 转 纳米
    pub fn angstrom_to_nm(value: f64) -> f64 {
        value * 0.1
    }
    
    /// 纳米 转 Angstrom
    pub fn nm_to_angstrom(value: f64) -> f64 {
        value * 10.0
    }
    
    /// Bohr 转 Angstrom
    pub fn bohr_to_angstrom(value: f64) -> f64 {
        value * 0.529177
    }
    
    /// Angstrom 转 Bohr
    pub fn angstrom_to_bohr(value: f64) -> f64 {
        value / 0.529177
    }
}

/// 能量单位转换
pub mod energy {
    /// Hartree 转 kcal/mol
    pub fn hartree_to_kcal_mol(value: f64) -> f64 {
        value * 627.509
    }
    
    /// kcal/mol 转 Hartree
    pub fn kcal_mol_to_hartree(value: f64) -> f64 {
        value / 627.509
    }
    
    /// eV 转 kcal/mol
    pub fn ev_to_kcal_mol(value: f64) -> f64 {
        value * 23.061
    }
    
    /// kcal/mol 转 kJ/mol
    pub fn kcal_mol_to_kj_mol(value: f64) -> f64 {
        value * 4.184
    }
    
    /// kJ/mol 转 kcal/mol
    pub fn kj_mol_to_kcal_mol(value: f64) -> f64 {
        value / 4.184
    }
}

/// 温度单位转换
pub mod temperature {
    /// Celsius 转 Kelvin
    pub fn celsius_to_kelvin(value: f64) -> f64 {
        value + 273.15
    }
    
    /// Kelvin 转 Celsius
    pub fn kelvin_to_celsius(value: f64) -> f64 {
        value - 273.15
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_length_conversion() {
        assert!((length::angstrom_to_nm(10.0) - 1.0).abs() < 1e-10);
        assert!((length::nm_to_angstrom(1.0) - 10.0).abs() < 1e-10);
    }
    
    #[test]
    fn test_energy_conversion() {
        let hartree = 1.0;
        let kcal = energy::hartree_to_kcal_mol(hartree);
        assert!((kcal - 627.509).abs() < 0.001);
    }
}
