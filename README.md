# DemStochLifeHisto
# Project Title

Simulating the difffusion approximations from "Effects of demographic stochasticity and life-history strategies on times and probabilities to fixation"

## Getting Started

Here we provide the C++ script and other relevant files used to simulate the bi-dimensional diffusion process given by equations 1a and 1b  in "Effects of demographic stochasticity and life-history strategies on times and probabilities to fixation". 

### Prerequisites

Three files are provided "main.cpp", "MersenneTwister.h" the header file containing random number generator required for the simulations and "param.txt" the text file required to input the parameters (an example of a parameter set is already in the param.txt file).


### Installing

In order to compile the script, a compiler such as GNU is required.
From the terminal a simple command can compile the script:

g++ -o ProgName main.cpp MersenneTwister.h -lm

## Running the simulations

Once compiled the command "./ProgName" (or "nohup ./ProgName >out 2>&1&" so as to rendre the simulation independent of the terminal) suffices to launch the program. Several parameter sets can be included in "param.txt", each begining with an "*" so as to be read by the script. They will be run in sequence.

## Parameters used
typesim: 0 Our model with genetic feedback, 1 Simulations without feed-back (population size remains stochastic however), and 2 Simulations with constant population size
mode: 0 writes only means for each of the n repetitions at the end of simulation, the last line representing the summary statistics over all n repetitions, 1 writes entire trajectory for 1 simulation only (n = 1);
N0: initial population size
X0: initial frequency of allele a
b: birth rate
d: death rate
c: competition
alpha: self-fertilisation rate
g: parameter gamma for scaling the speed of birth and death rates. In order to obtain a model as close as possible to the Wright-Fisher diffusion, this parameter is set to 0.5.
s: selection coefficient (sigma in the paper)
h: dominance
dt: Time step, set to 10e-04 for large populaiton and growth rates, and 10e-05 otherwise
n: number of repeats for each parameter set
id: identity number defined by the used and added at the end of the output file so as to avoid over-writing simulation outputs with the same parameter values