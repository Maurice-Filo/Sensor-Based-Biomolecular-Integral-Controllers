%% Clear Workspace
close all
clear
clc
disp(['Running: "', mfilename, '.m"']);

%% Configure Working Environment
ProjectPath = '/cluster/home/mfilo/Fine Exhaustive Code';
DataPath = fullfile(ProjectPath, 'Data');
addpath(genpath(ProjectPath));
euler = parcluster('local');
pool = parpool(euler, 16);

% ProjectPath = '/Users/mfilo/Library/CloudStorage/GoogleDrive-maurice.g.filo@gmail.com/My Drive/Research at ETH/My Papers and Presentations/Bacterial PI/Matlab/Noise/Fine Exhaustive Code';
% DataPath = fullfile(ProjectPath, 'Data');
% addpath(genpath(ProjectPath));

%% Network Construction
StoichiometryMatrix = StoichiometryMatrix_FP_BirthDeath();
PropensityFunction = @PropensityFunction_FP_BirthDeath;
Parameters = Parameters_FP_BirthDeath();
Parameters.theta = 10;
Output_Index = 1;

%% Simulation Settings
IC = [0; 0];
TimeSpan_Long = [0, 1e6];
N_Trajectories_Few = 1;
TransientPortion = 0.01;

%% Swept Parameters
N_delta = 64;
delta_vector = logspace(-3, 4, N_delta);

%% Simulations
SimTime = tic;
StationaryMean = zeros(N_delta, 1);
StationaryVariance = zeros(N_delta, 1);
disp(['Progress: (i) = (:) out of (', num2str(N_delta), ')']);
parfor i = 1 : N_delta
    localParameters = Parameters;
    localParameters.delta = delta_vector(i);
    [StationaryMeanX, StationaryVarianceX, ~] = ComputeStatistics_Stationary(PropensityFunction, StoichiometryMatrix, localParameters, IC, TimeSpan_Long, TransientPortion, N_Trajectories_Few);
    StationaryMean(i) = StationaryMeanX(Output_Index);
    StationaryVariance(i) = StationaryVarianceX(Output_Index);
end
SimTime = toc(SimTime);
disp(['Successfully completed the run of: "', mfilename, '.m". Simulation Time is: ', num2str(SimTime/3600), ' hours.']);

%% Save data
Save_Flag = 0;
FileName = fullfile(DataPath, 'FP_BirthDeath_FixedHightheta');
if Save_Flag == 1
    clear euler pool
    save(FileName);
end