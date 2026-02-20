# WOTC Code Repository

This repository accompanies the thesis and serves as a guide for reproducing simulations and generating the figures used in the text. You can find the source repository at:
https://github.com/lenhoy/WOTC

## Dependencies

* MATLAB and Simulink
* MSS Toolbox for Simulink

## Getting Started

1. Clone or download the repository.
2. Open MATLAB and set the current folder to the repository root.
3. Open the project file (`WOTC.prj`) to set up paths.
4. Open the main Simulink model (`WOTC_simulation.slx`).
5. Select a scenario using the `Signal Editor` block within the `Guidance` subsystem. The folder `Scenarios` contains the relevant files.
6. Select a controller. Inside the `Controller` subsystem, a clickable switch can be used to toggle between the `PID` and `WOTC` controllers.
7. Settings.
8. Signals of interest are logged and can be viewed using the `Data Inspector`.

## Repository Structure

The following is a quick description of relevant folders used by the simulation or for saving data and for plotting purposes.

### Callbacks
Before each variables are cleared and `run('models/osvData.m');` is used to reset the workspace. After each run the stop callback saves the output data. 

### Output Data
The simulation data outputs are stored in the `output` data folder. A post-run callback automatically saves timestamped data for each run (`.mat` files).

### Scenarios
The `scenarios` folder contains the following:
* Scenario files used by the Signal Editor.
* Script files to generate trajectory data.

## Viewing and Plotting Results

**Recommended: Simulink Data Inspector**
For interactive analysis and quick comparisons. Most signals are already logged and are available to view after a run.

**Python Plots**
The repository also contains a `plots` folder with Python scripts that reproduce the same plots as shown in the thesis.

1. Run a simulation and locate the desired timestamped `.mat` file in the output folder.
2. Convert the output file to Python using:
   `convert_for_python.m`
   *(Note: This sometimes needs to be added to your MATLAB path.)*
3. In `plots/python/scripts` a `readme` outlines the available plots.

---
*Note: AI tools were used in the creation and formatting of this README.*
