# STL_STT: Approximation-free Control for Signal Temporal Logic Specifications using Spatiotemporal Tubes

This repository provides the implementation for synthesizing spatiotemporal tubes (STTs) and designing controllers to satisfy Signal Temporal Logic (STL) specifications. 
This has been presented in the paper "Approximation-free Control for Signal Temporal Logic Specifications using Spatiotemporal Tubes"

It includes two case studies:

- **Omnidirectional Robot**
- **Rotating Rigid Spacecraft**

## Repository Structure

### Tube Generation (Python)
- `OmnidirectionalRobot.py`: Generates STTs for the omnidirectional robot example.
- `RotatingRigidSpacecraft.py`: Generates STTs for the rotating rigid spacecraft example.

### Control and Visualization (MATLAB)
- `Omnirobot.m`: Synthesizes the closed-form controller and simulates the omnidirectional robot trajectory using the generated tube.
- `space.m`: Synthesizes the closed-form controller and simulates the spacecraft trajectory using the generated tube.
- `plotcube.m`: Utility script for visualizing 3D tube structures.

### Data and Results
- `Robot.csv`, `Robot1.csv`: Tube and trajectory data for the omnidirectional robot.
- `Spacecraft.csv`, `Spacecraft1.csv`, `Spacecraft2.csv`: Tube and trajectory data for the spacecraft.
- `Omni2.fig`, `Spacecraft_1.fig`: MATLAB figure files for visualizations.
- `Omni2.png`, `Spacecraft_1.png`: Plots of the robot and spacecraft trajectories and STTs.

## Requirements

- Python 3.7+ with `z3-solver` installed (`pip install z3-solver`)
- MATLAB R2024B or later (for running `.m` files)

## How to Run

1. **Generate STTs**:
   - Run `OmnidirectionalRobot.py` or `RotatingRigidSpacecraft.py` in Python to generate STTs.

2. **Run Controller and Simulation**:
   - Open MATLAB.
   - Run `Omnirobot.m` for the robot example or `space.m` for the spacecraft.
   - This will simulate the system within the tube and visualize the results.

## Citation

If you use this code in your research, please cite the following paper:

```bibtex
@article{STT_STL_arxiv,
  title={Approximation-free Control for Signal Temporal Logic Specifications using Spatiotemporal Tubes},
  author={Das, Ratnangshu and Choudhury, Subhodeep and Jagtap, Pushpak},
  journal={arXiv preprint arXiv:2505.05323},
  year={2024}
}