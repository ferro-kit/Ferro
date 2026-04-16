//! 轨迹数据结构

use crate::Molecule;
use serde::{Serialize, Deserialize};

/// 分子动力学轨迹
/// 
/// 存储多个时间步的分子结构
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Trajectory {
    /// 每个时间步的分子结构
    pub frames: Vec<Molecule>,
    /// 时间步长 (单位: ps)
    pub timestep: Option<f64>,
    /// 总时间 (单位: ps)
    pub total_time: Option<f64>,
}

impl Trajectory {
    /// 创建空轨迹
    pub fn new() -> Self {
        Self {
            frames: Vec::new(),
            timestep: None,
            total_time: None,
        }
    }
    
    /// 添加一帧
    pub fn add_frame(&mut self, frame: Molecule) {
        self.frames.push(frame);
    }
    
    /// 获取帧数
    pub fn frame_count(&self) -> usize {
        self.frames.len()
    }
    
    /// 获取指定帧
    pub fn get_frame(&self, index: usize) -> Option<&Molecule> {
        self.frames.get(index)
    }
    
    /// 获取第一帧
    pub fn first_frame(&self) -> Option<&Molecule> {
        self.frames.first()
    }
    
    /// 获取最后一帧
    pub fn last_frame(&self) -> Option<&Molecule> {
        self.frames.last()
    }
    
    /// 设置时间步长
    pub fn set_timestep(&mut self, dt: f64) {
        self.timestep = Some(dt);
        if !self.frames.is_empty() {
            self.total_time = Some(dt * (self.frames.len() - 1) as f64);
        }
    }
    
    /// 获取指定帧的时间
    pub fn get_time(&self, frame_index: usize) -> Option<f64> {
        self.timestep.map(|dt| dt * frame_index as f64)
    }
}

impl Default for Trajectory {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Atom;
    use nalgebra::Vector3;
    
    #[test]
    fn test_trajectory_creation() {
        let traj = Trajectory::new();
        assert_eq!(traj.frame_count(), 0);
    }
    
    #[test]
    fn test_add_frames() {
        let mut traj = Trajectory::new();
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new("C", Vector3::new(0.0, 0.0, 0.0)));
        
        traj.add_frame(mol.clone());
        traj.add_frame(mol);
        
        assert_eq!(traj.frame_count(), 2);
    }
}
