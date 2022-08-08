# Rigorous data-driven computation of spectral properties of Koopman operators for dynamical systems

## Abstract
When dealing with dynamical systems with discrete time steps, nonlinear and unknown dynamics 
pose challenges to the standard Poincaré’s approach. Koopman operators are infinite-dimensional
operators that linearize the dynamics by lifting it to the observables state space. We review spectral properties of Koopman operators, focusing in particular on the case of linear dynamics, as well
as numerical methods for computing spectral information of the Koopman operator from snapshot
data. We also examine the problem of spectral pollution, recalling Residual Dynamic Mode Decomposition (ResDMD) and proving a proposition that allows reducing ResDMD to a standard
matrix pseudospectrum approximation problem. Numerical examples from [1] are reproduced. Finally, a two step procedure for extracting a dictionary of observables directly from snapshot data
and a general-purpose kernel is introduced.

### Main Reference
[1] Matthew J. Colbrook and Alex Townsend, [Rigorous data-driven computation of spectral properties of Koopman operators for dynamical systems](https://doi.org/10.48550/arXiv.2111.14889), Nov. 2021.

## Requirements
MATLAB R2019b or later versions.

## Repository Description
* `code` - Folder containing code and auxiliary files. The file `run.m` reproduces all the figures in the report.
* `presentation` - Folder containing the slides for presenting the project
* `report` - Folder containing the report

## Authors
- Author: Ivan Bioli
- Professor: D. Kressner
- Supervisor: Alice Cortinovis