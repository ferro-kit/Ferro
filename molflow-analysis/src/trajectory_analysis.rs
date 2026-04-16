//! 轨迹分析工具

use molflow_core::Trajectory;
use nalgebra::Vector3;
use rayon::prelude::*;

/// 计算均方位移 (MSD)
pub fn mean_squared_displacement(traj: &Trajectory) -> Vec<f64> {
    if traj.frames.is_empty() {
        return vec![];
    }
    
    let n_frames = traj.frame_count();
    let n_atoms = traj.first_frame().map(|f| f.atom_count()).unwrap_or(0);
    
    let mut msd = vec![0.0; n_frames];
    
    for dt in 0..n_frames {
        let mut sum = 0.0;
        let mut count = 0;
        
        for t0 in 0..(n_frames - dt) {
            let frame0 = &traj.frames[t0];
            let frame1 = &traj.frames[t0 + dt];
            
            for i in 0..n_atoms {
                if let (Some(atom0), Some(atom1)) = (frame0.atoms.get(i), frame1.atoms.get(i)) {
                    let displacement = atom1.position - atom0.position;
                    sum += displacement.norm_squared();
                    count += 1;
                }
            }
        }
        
        if count > 0 {
            msd[dt] = sum / count as f64;
        }
    }
    
    msd
}

/// 计算均方根偏差 (RMSD) 相对于第一帧
pub fn rmsd_to_first(traj: &Trajectory) -> Vec<f64> {
    if traj.frames.is_empty() {
        return vec![];
    }
    
    let reference = &traj.frames[0];
    let n_atoms = reference.atom_count();
    
    traj.frames.par_iter()
        .map(|frame| {
            let mut sum = 0.0;
            for i in 0..n_atoms {
                if let (Some(ref_atom), Some(frame_atom)) = 
                    (reference.atoms.get(i), frame.atoms.get(i)) {
                    let diff = frame_atom.position - ref_atom.position;
                    sum += diff.norm_squared();
                }
            }
            (sum / n_atoms as f64).sqrt()
        })
        .collect()
}

/// 计算分子质心随时间的轨迹
pub fn center_of_mass_trajectory(traj: &Trajectory) -> Vec<Vector3<f64>> {
    traj.frames.iter()
        .map(|frame| frame.center_of_mass())
        .collect()
}

/// 计算回转半径随时间的变化
pub fn radius_of_gyration_trajectory(traj: &Trajectory) -> Vec<f64> {
    traj.frames.par_iter()
        .map(|frame| {
            let com = frame.center_of_mass();
            let total_mass = frame.total_mass();
            
            if total_mass == 0.0 {
                return 0.0;
            }
            
            let sum: f64 = frame.atoms.iter()
                .map(|atom| {
                    let r = (atom.position - com).norm();
                    atom.mass * r * r
                })
                .sum();
            
            (sum / total_mass).sqrt()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use molflow_core::{Molecule, Atom};
    
    #[test]
    fn test_trajectory_analysis() {
        let mut traj = Trajectory::new();
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new("C", Vector3::new(0.0, 0.0, 0.0)));
        
        traj.add_frame(mol.clone());
        traj.add_frame(mol);
        
        let rmsd = rmsd_to_first(&traj);
        assert_eq!(rmsd.len(), 2);
    }
}
