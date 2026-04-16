//! ChemTools IO - 文件读写模块
//! 
//! 支持多种分子结构文件格式的读取和写入

pub mod readers;
pub mod writers;

// 重新导出常用函数，避免命名冲突
pub use readers::{read_xyz, read_pdb};
pub use writers::{write_xyz, write_pdb};

#[cfg(test)]
mod tests {
    #[test]
    fn test_io_module() {
        // 基本测试确保模块可以编译
        assert!(true);
    }
}
