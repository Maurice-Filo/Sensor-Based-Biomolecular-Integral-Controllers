%% Clear Workspace
% close all
clear
clc
disp(['Running: "', mfilename, '.m"']);

%% Configure Working Environment
% ProjectPath = '/cluster/home/mfilo/Code';
% DataPath = fullfile(ProjectPath, 'Data');
% addpath(genpath(ProjectPath));
% euler = parcluster('local');
% pool = parpool(euler, 48);

ProjectPath = '/Users/mfilo/Library/CloudStorage/GoogleDrive-maurice.g.filo@gmail.com/My Drive/Research at ETH/My Papers and Presentations/Bacterial PI/Matlab/Noise/Code';
DataPath = fullfile(ProjectPath, 'Data');
addpath(genpath(ProjectPath));

%% Network Construction
StoichiometryMatrix = StoichiometryMatrix_sAIF_BirthDeath();
PropensityFunction = @PropensityFunction_sAIF_BirthDeath;
Parameters = Parameters_sAIF_BirthDeath();
Output_Index = 1;

%% Parameters
Parameters.delta = 0.1;
Parameters.theta = 10;
Parameters.mu = 0;
Parameters.eta = 0;

%% Simulation Settings
IC = [0; 0; 0];
TimeSpan_Long = [0, 1e6];
N_Trajectories_Few = 1;
TransientPortion = 0.01;

TimeSpan = [0, 300];
N_Trajectories = 1e2;
N_Grid = 100;

%% Simulation using Time Averaging
SimTime1 = tic;
[StationaryMeanX, StationaryVarianceX, ~] = ComputeStatistics_Stationary(PropensityFunction, StoichiometryMatrix, Parameters, IC, TimeSpan_Long, TransientPortion, N_Trajectories_Few);
StationaryMean = StationaryMeanX(Output_Index);
StationaryVariance = StationaryVarianceX(Output_Index);
StationaryCV = sqrt(StationaryVariance) / StationaryMean;
SimTime1 = toc(SimTime1);
disp(['Successfully completed the run of: "', mfilename, '.m" using time averaging. Simulation Time is: ', num2str(SimTime1), ' seconds.']);

%% Simulation using Ensemble Averaging
SimTime2 = tic;
[T, X, ~] = GenerateTrajectories_Grid(StoichiometryMatrix, PropensityFunction, Parameters, IC, TimeSpan, N_Trajectories, N_Grid);
[MeanX, VarX] = ComputeStatistics_Grid(T, X);
Mean = MeanX(Output_Index,:); 
Variance = VarX(Output_Index,:);
CV = sqrt(Variance) ./ Mean;
SimTime2 = toc(SimTime2);
disp(['Successfully completed the run of: "', mfilename, '.m" using ensemble averaging. Simulation Time is: ', num2str(SimTime2), ' seconds.']);

%% Print Results
disp(Parameters)
disp(['Time averaging: Mean = ', num2str(StationaryMean), '; CV = ', num2str(StationaryCV)]);
disp(['Ensemble averaging: Mean = ', num2str(Mean(end)), '; CV = ', num2str(CV(end))]);

%% Plot Results
figure();
subplot(1,2,1);
plot(T{1}, Mean, 'Color', 'k');
hold on;
plot(T{1}(end), StationaryMean, 'Marker', 'o', 'Color', 'r');
subplot(1,2,2);
plot(T{1}, CV, 'Color', 'k');
hold on;
plot(T{1}(end), StationaryCV, 'Marker', 'o', 'Color', 'r');