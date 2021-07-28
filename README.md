# Bayesian nonparametric inference for heterogeneously-mixing infectious disease models
---
The repositroy contains code for simulating outbreaks of heterogneously mixing infectious diseases and fitting models to outbreak data. In the model, the infection rate from individual _i_ to _j_ is _beta{ij}_. The methodology is Bayesian nonparametric meaning the functional form of _beta{ij}_ does not need to be specified. The form is inferred using a Gaussian process. The model allows for the population to be split into a number of types, e.g. cattle farms and pig farms, and the infection rate to depend on the type of infectious farm. In cases where the times individuals are infected are unknown, these can be inferred using the code provided. 

---

## Overview

### C code
The code avaliable in C provides a multitude of modelling options and can be tailored to a wide variety of circumstances. It is effiecient and can supports simulation for very large outbreaks. It requires compiltation using `gcc`. It can be used with `openmp` to speed up the computations, but this requires multi-core computation. Simulation and inference is made easier through a bash script where the user can specify the variables relating to the epidemic. The script then simulates, performs inference for an epidemic, and plots the inferred results. 

### R Code
The R code provided is the most user friendly code in this repository and provides a good starting place for the user. Outbreaks can be simulated using the R code and the inference code provides a implementation of the single type GP model. Due to the limitations of R only small and medium size outbreaks are supported (N < 2000). For larger outbreaks or where infection times are unknown, the C code should be used. 

---

## Models

### Multi-Output covariance model
In this model the infection rate functions for each type are correlated. The allows the relationship between each pair of types to be summarised through one parameter. The model should be run using the bash script provided. It simulates locations for N individuals and assigns them each one a P types, where N and P are specified by the user. The functional form for each type is exponential, although this can be easily changed. The inference code performs inferences for this data or any user provided data set which contains coordinates, type and a removal time for each individual. To use the independent model, the correlation parameter should be set to 0. 

### Discrepancy model
This models selects one type as a baseline and then computes the difference between the infection rate between each type and the baseline type. This makes it easier to interpret and visualise results. The type which has label 1 should be taken as the baseline. The inference code performs inferences for this data or any user provided data set which contains coordinates, type and a removal time for each individual.

### Single type model
The single type model assumes all types are equal and the infection rate only depends on a function of the pair of individuals, e.g. pair-wise distance. This is implemented in R, so provides a user friendly introduction to the inference algorithm. 

