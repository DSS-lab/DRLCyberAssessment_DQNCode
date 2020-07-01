# DeSAW Code

Repository containing supplementary test cases and code for the paper "DeSAW: Deep Reinforcement Learning for Security Assessment of Wind Integrated Power Systems" by Xiaorui Liu, Juan Ospina, Darbi Bauer, and Charalambos Konstantinou.



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

## Running the code

All necessary functions are provided by the file test.jl. The cases folder includes Matpower cases which could be directly used as test cases in Julia. All the required files for POMDP.jl are in the folder "extra".





