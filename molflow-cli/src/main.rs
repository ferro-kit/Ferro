//! ChemTools 命令行工具
//! 
//! 提供便捷的命令行接口来使用 ChemTools 库的各种功能

use clap::{Parser, Subcommand};
use anyhow::Result;

#[derive(Parser)]
#[command(name = "molflow")]
#[command(about = "计算化学工具集", long_about = None)]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// 文件格式转换
    /// 示例: molflow convert -i input.xyz -o output.pdb
    Convert {
        /// 输入文件
        #[arg(short, long)]
        input: String,
        
        /// 输出文件
        #[arg(short, long)]
        output: String,
    },
    
    /// 分析轨迹
    /// 示例: molflow analyze -t traj.xyz -p rmsd
    Analyze {
        /// 轨迹文件
        #[arg(short, long)]
        traj: String,
        
        /// 分析类型 (rmsd, msd, rg)
        #[arg(short, long)]
        property: String,
    },
    
    /// 创建计算任务
    /// 示例: molflow job -i input.xyz -s gaussian -m B3LYP
    Job {
        /// 输入结构文件
        #[arg(short, long)]
        input: String,
        
        /// 软件类型 (gaussian, orca, gromacs)
        #[arg(short, long)]
        software: String,
        
        /// 计算方法
        #[arg(short, long)]
        method: Option<String>,
        
        /// 输出文件
        #[arg(short, long)]
        output: Option<String>,
    },
    
    /// 显示分子信息
    /// 示例: molflow info -i molecule.xyz
    Info {
        /// 输入文件
        #[arg(short, long)]
        input: String,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    
    match cli.command {
        Commands::Convert { input, output } => {
            println!("转换 {} -> {}", input, output);
            // 调用 molflow-io 进行格式转换
            convert_file(&input, &output)?;
        }
        
        Commands::Analyze { traj, property } => {
            println!("分析轨迹 {} 的性质: {}", traj, property);
            // 调用 molflow-analysis 进行分析
            analyze_trajectory(&traj, &property)?;
        }
        
        Commands::Job { input, software, method, output } => {
            println!("从 {} 创建 {} 任务", input, software);
            // 调用 molflow-workflow 创建任务
            create_job(&input, &software, method.as_deref(), output.as_deref())?;
        }
        
        Commands::Info { input } => {
            println!("显示 {} 的信息", input);
            // 调用 molflow-io 和 molflow-core 显示信息
            show_info(&input)?;
        }
    }
    
    Ok(())
}

fn convert_file(input: &str, output: &str) -> Result<()> {
    use molflow_io::{read_xyz, write_pdb};
    
    // 简单示例：读取 XYZ，写入 PDB
    let mol = read_xyz(input)?;
    write_pdb(&mol, output)?;
    println!("转换完成!");
    
    Ok(())
}

fn analyze_trajectory(traj: &str, property: &str) -> Result<()> {
    // TODO: 实现轨迹分析
    println!("轨迹分析功能待实现: {} {}", traj, property);
    Ok(())
}

fn create_job(input: &str, software: &str, method: Option<&str>, output: Option<&str>) -> Result<()> {
    use molflow_io::read_xyz;
    use molflow_workflow::GaussianJobBuilder;
    
    let traj = read_xyz(input)?;
    let frame = traj.frames.into_iter().next()
        .ok_or_else(|| anyhow::anyhow!("文件中没有帧"))?;

    match software {
        "gaussian" => {
            let mut builder = GaussianJobBuilder::new(frame);
            if let Some(m) = method {
                builder.method = m.to_string();
            }
            let job_file = builder.build()?;
            
            let output_path = output.unwrap_or("job.gjf");
            std::fs::write(output_path, job_file)?;
            println!("已创建 Gaussian 任务文件: {}", output_path);
        }
        _ => {
            println!("暂不支持的软件: {}", software);
        }
    }
    
    Ok(())
}

fn show_info(input: &str) -> Result<()> {
    use molflow_io::read_xyz;
    
    let traj = read_xyz(input)?;
    let frame = traj.frames.first()
        .ok_or_else(|| anyhow::anyhow!("文件中没有帧"))?;
    println!("分子信息:");
    println!("  帧数:   {}", traj.frames.len());
    println!("  原子数: {}", frame.atoms.len());
    
    Ok(())
}
