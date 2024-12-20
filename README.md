# Simulation of Discrete Shell behavior on a Paper Mesh

This project implements the paper Discrete Shells from Grinspun et al. where a precise physical behavior of thin shells got simulated.
We use libigl for the representation of the mesh in a gui and several other libigl functions for modifying the state of the mesh.
Additionally automatic differentiation tools from TinyAD are used for simplifying several derivation calculations.

## Project Structure
```
thin_shells/
├── CMakeLists.txt
├── src/
│   ├── constants.h              # Header file defining project constants
│   ├── DiscreteShell.cpp        # Implementation of the Discrete Shell class
│   ├── DiscreteShell.h          # Header file for the Discrete Shell class
│   ├── main.cpp                 # Main program entry point
│   ├── Mesh.cpp                 # Implementation of the Mesh class
│   ├── Mesh.h                   # Header file for the Mesh class
│   ├── Solver.cpp               # Implementation of the Solver class
│   ├── Solver.h                 # Header file for the Solver class
│   └── utils.h                  # Utility functions for common tasks
├── tinyAD/                 # Contain external libraries for automatic differentiation
├── stb/                    # Contain external libraries for .PNG related code
├── autodiff/               # Contain external libraries for automatic differentiation
├── build/
│   └── ...
├── docs/               # Contain used papers or a summary of the paper used
└── README.md
```
## Installation Windows
1. Clone this repository:
```
git clone <repository_url>
cd <repository_name>
```
2. Build the project using cmake GUI (working version 3.31.2):
3. Configure and Generate 
4. Open thinshells.sln with Visual Studio 2017 or higher
5. In Visual Studio open configure Startup Projects, and define thinshells.sln file as single startup project.

## Installation Linux
```
mkdir build
cd build
cmake ..
make 
```

## Usage Windows
To run the program, use the following commands:
```
1. Use Visual Studio Interface to start up the project in Release mode.
2. press " " to start the simulation.
```
To run the program with a specific mesh:
```
1. in main.cpp line 23 change the string filename to a certain "<filename>.off" file in the data folder.
ATTENTION!: read important information first
2. press " " to start the simulation.
```

### Arguments
First argument also defines the "<filename>.off" which is to be run, otherwise it just takes a simple flat rectangle mesh as starting point.


## Important information

The simulation does not work for meshes which are not flat when reading the file. 
I.e. twisted.off and rectangle_folded.off are examples which do not work.

