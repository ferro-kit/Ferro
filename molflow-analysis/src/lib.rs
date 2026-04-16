//! ChemTools Analysis - 后处理分析模块
//! 
//! 提供轨迹分析、几何计算、性质计算等功能

pub mod geometry;
pub mod trajectory_analysis;
pub mod properties;

pub use geometry::*;
pub use trajectory_analysis::*;
pub use properties::*;

#[cfg(test)]
mod tests {
    #[test]
    fn test_analysis_module() {
        assert!(true);
    }
}
