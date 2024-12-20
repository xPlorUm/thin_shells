# Simulation of Discrete Shell behavior on a Paper Mesh

This project implements the paper Discrete Shells from Grinspun et al. where a precise physical behavior of thin shells got simulated.
We use libigl for the representation of the mesh in a gui and several other libigl functions for modifying the state of the mesh.
Additionally automatic differentiation tools from TinyAD are used for simplifying several derivation calculations.

## Project Structure
```
src/
│
├── datasets/
│   └── load_datasets.py          # Dataset loading and task splitting logic
│
├── models/
│   └── resnet.py                 # Define ResNet or other models
│
├── replay_buffer/
│   └── fifo_replay_buffer.py     # ReplayBuffer class (shared across all strategies)
│
├── training_strategies/
│   ├── base_strategy.py             # Base training class
│   ├── goldilock_strategy.py        # Goldilock-based training
│   ├── goldilock_v2_strategy.py     # Fast Goldilock-based training using predicted learning-speeds 
│   ├── weighted_buffer_strategy.py  # Weighted Replaybuffer training
│   └── replaybuffer_strategy.py     # Replaybuffer based training
│
├── evaluation_methods/
│   └── basic_plots.py            # Evaluation scripts for testing performance
│
├── results/
│   └── logs/                     # Training logs and metrics
│   └── plots/                    # Visualization of results
│
├── requirement.txt               # Main entry point for the project
├── main.py                       # Main entry point for the project
└── README.md                     # Project description and usage (This file!)
```
## Features
- Replay Buffer: Implements a replay buffer to store and sample data for training to help with catastrophic forgetting in continual learning scenarios.
- Multiple Training Strategies:
    - BaseStrategy: A basic strategy for training the model without using a replay buffer.
    - ReplayBufferStrategy: A strategy that uses the replay buffer to store past experiences and alleviate forgetting.
    - GoldilockStrategy: A strategy based on the learning speed of different tasks to determine which data should be added to the buffer.
    - WeightedStrategy: A strategy based on the learning speed. slower learning speed examples have higher sampling weight
- Models: The project includes the ResNet18 and ResNet34 architectures for the model.
- Logging: All results and training logs are saved in the results/log/ directory.
- Evaluation: The model is evaluated based on test loss and accuracy after each task and training loss for each epoch.

## Installation
1. Clone this repository:
```
git clone <repository_url>
cd <repository_name>
```
2. Build the project using cmake:
```
mkdir build
cd build
cmake ..
make 
```

## Usage
To run the program, use the following commands:
```

```
To run the program with a specific mesh:
```
```
### Arguments
- `--dataset`: The dataset to use. Options: CIFAR-10, CIFAR-100, TinyImagenet.



