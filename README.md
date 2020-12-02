# Code supplement for "DRL for Cybersecurity Assessment of Wind Integrated Power Systems"


This repository contains supplementary test cases and code for the paper "Deep Reinforcement Learning for Cybersecurity Assessment of Wind Integrated Power Systems", by XiaoRui Liu, Juan Ospina, and Charalambos Konstantinou. 

If you find this code useful, please consider citing our paper:

@ARTICLE{9261465,
  author={X. {Liu} and J. {Ospina} and C. {Konstantinou}},
  journal={IEEE Access}, 
  title={Deep Reinforcement Learning for Cybersecurity Assessment of Wind Integrated Power Systems}, 
  year={2020},
  volume={8},
  number={},
  pages={208378-208394},
  doi={10.1109/ACCESS.2020.3038769}}

## Julia Installation

You can download Julia from [Link](https://julialang.org/).


## Getting started with Julia

The code and packages are in the "Julia" programming language. 

Before running the code, make sure all packages are installed. The packages can be added in a Julia environment directly via any source code editor which supports "Julia". Here we use "Atom". 


```python
import packages:

julia> Pkg.add("POMDPSimulators")
julia> Pkg.add("RLInterface")
julia> Pkg.add("POMDPModelTools")

....

```

## Version of Packages 

```bash
Atom v0.11.2
BSON v0.2.4
DeepQLearning v0.3.0
FileIO v1.0.7
Flux v0.9.0
Ipopt v0.6.0
Json v0.21.0
Juno v0.7.2
MATLAB v0.7.3
POMDPModelTools v0.2.0
POMDPModels v0.4.0
POMDPSimulators v0.3.1
POMDPS v0.8.1
PowerModels v0.13.0
```



## Running the code

All necessary functions are provided by the file test.jl. The cases folder includes Matpower cases which could be directly used as test cases in Julia. All the required files for POMDP.jl are in the folder "extra".


