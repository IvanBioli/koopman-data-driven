clear all
addpath(genpath(pwd))
% NOTE: Figures 3.2.b and 3.3.b require a long runtime. In order to load
% the results from the provided workspaces, the flag load is set to true.
% Should you want to run the code entirely, please set this flag to false.
loading = true;
% 2.3: Comparison of the convergence for DMD, SVD-based DMD and Arnoldi
DMD_vs_Arnoldi;
% 3.3.1: Gauss iterated map
Gauss_iterated_map;
% 4: Kernalized EDMD
Kernelized_Gauss_iterated_map;
% 3.3.2: Nonlinear pendulum
N = 20;
Nonlinear_pendulum;
N = 100;
Nonlinear_pendulum;