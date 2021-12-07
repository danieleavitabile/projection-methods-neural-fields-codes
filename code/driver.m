%% Cleaning
clear all, close all, clc;

%% Fetch current folder
rootPath = pwd;

%% Results folder
resultsPath = [rootPath '/../results'];

%% Running tests
cd([rootPath '/FECollocation/Trapezium']); RunTest(rootPath,resultsPath); cd(rootPath);
cd([rootPath '/SpectralCollocation/Trapezium']); RunTest(rootPath,resultsPath); cd(rootPath);
cd([rootPath '/SpectralCollocation/ClenshawCurtis']); RunTest(rootPath,resultsPath); cd(rootPath);
cd([rootPath '/FEGalerkin/Gauss']); RunTest(rootPath,resultsPath); cd(rootPath);
cd([rootPath '/SpectralGalerkin/Trapezium']); RunTest(rootPath,resultsPath); cd(rootPath);
