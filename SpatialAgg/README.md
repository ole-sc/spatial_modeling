## Usage Instructions 
- To run the code, you have to install julia (obviously).
- Navigate to the directory of the module `cd path/to/SpatialAgg`.
- type `julia` to start the REPL.
- hit `]` to open the package manager.
- Activate the module `activate .` 
- Load dependencies (still in the package manager) `instantiate`and wait for the installation to finish (this might take a while).
- exit the package manager with `ctrl+c`.
- Open the file "tests/runexperiments" in a text-editor and uncomment the lines to run certain pieces of code (there are three real-time plots that are recommended and then all other simulations I ran to generate the figures for the report).
- In the julia REPL, type `include("tests/runexperiments.jl")` to run the simulation (It will probably take a moment before anything happens, because the first run in julia is always slow).

## Files
- `src` contains the source files for the simulation
- `tests` contains files for running experiments, testing codes and making plots.
- `plots` has the heatmaps for N=200 and N=5000 that were not in the report.