# Implicit Finite Difference Scheme for Acoustic Wave Propagation 🌊🔬

A numerical simulation project using Python and MATLAB to model acoustic wave propagation using implicit finite difference and central difference schemes.

## Table of Contents

* [Project Overview](#project-overview)
* [Dataset / Inputs](#dataset--inputs)
* [Features](#features)
* [Installation](#installation)
* [Usage](#usage)
* [Implementation](#implementation)
* [Results](#results)
* [Future Work](#future-work)
* [License](#license)

## Project Overview

This project models acoustic wave propagation using implicit finite difference and central difference schemes. By simulating waves at different propagation speeds using the Richter wavelet, the project provides insights into numerical stability, wave behavior, and computational methods for acoustic modeling.

Key components include:

* Implementation of implicit and central finite difference schemes.
* Simulation of wave propagation at four different speeds.
* Visualization of wave propagation using Python and MATLAB.

## Dataset / Inputs

* **Input Wavelet:** Richter wavelet used for generating initial acoustic signals.
* **Wave Propagation Speeds:** Simulations conducted for four different speeds to analyze numerical behavior.
* **Grid Parameters:** Spatial and temporal discretization parameters for the finite difference schemes.

## Features

* **Implicit Finite Difference Scheme:** Ensures numerical stability for simulating wave propagation.
* **Central Difference Scheme:** Provides baseline comparison for wave behavior.
* **Visualization:** Python scripts generate plots of wave propagation at different speeds.
* **MATLAB Integration:** MATLAB scripts used for additional visualization and analysis.

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/acoustic-wave-fd.git
cd acoustic-wave-fd

# Install Python dependencies
pip install -r requirements.txt
```

## Usage

1. Run Python scripts for visualization:

```bash
python simulate_wave.py
```

2. Open MATLAB scripts for further analysis:

```matlab
wave_propagation_simulation.m
```

3. Adjust parameters like wave speed, time steps, and spatial resolution to observe different behaviors.

## Implementation

* **Finite Difference Schemes:** Implemented implicit and central schemes for numerical simulation.
* **Wavelet Input:** Richter wavelet used as the source signal.
* **Visualization:** Python used for plotting temporal and spatial wave propagation graphs; MATLAB used for additional numerical analysis and verification.

## Results

* Successfully simulated acoustic wave propagation for four different speeds.
* Plotted wave profiles using both implicit and central schemes.
* Verified stability and accuracy of the implicit scheme compared to the central scheme.

## Future Work

* Extend simulation to 2D and 3D wave propagation models.
* Implement GPU acceleration for faster computation.
* Explore advanced wavelets and boundary conditions for realistic modeling.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
