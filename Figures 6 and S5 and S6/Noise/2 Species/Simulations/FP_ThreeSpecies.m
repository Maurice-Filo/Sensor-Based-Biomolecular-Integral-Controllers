%% Clear Workspace
close all
clear
clc
disp(['Running: "', mfilename, '.m"']);

%% Configure Working Environment
ProjectPath = '/cluster/home/mfilo/Fine Exhaustive Code Three Species';
DataPath = fullfile(ProjectPath, 'Data');
addpath(genpath(ProjectPath));
euler = parcluster('local');
pool = parpool(euler, 16);

% ProjectPath = '/Users/mfilo/Library/CloudStorage/GoogleDrive-maurice.g.filo@gmail.com/My Drive/Research at ETH/My Papers and Presentations/Bacterial PI/Matlab/Noise/Fine Exhaustive Code Three Species';
% DataPath = fullfile(ProjectPath, 'Data');
% addpath(genpath(ProjectPath));

%% Network Construction
StoichiometryMatrix = StoichiometryMatrix_FP_ThreeSpecies();
PropensityFunction = @PropensityFunction_FP_ThreeSpecies;
Parameters = Parameters_FP_ThreeSpecies();
Output_Index = 3;

%% Simulation Settings
IC = [0; 0; 0; 0];
TimeSpan_Long = [0, 1e6];
N_Trajectories_Few = 1;
TransientPortion = 0.01;

%% Swept Parameters
N_delta = 16;
N_theta = 32;
delta_vector = logspace(-3, -1, N_delta);
theta_vector = logspace(-5, 1, N_theta);

%% Simulations
SimTime = tic;
StationaryMean = zeros(N_delta, N_theta);
StationaryVariance = zeros(N_delta, N_theta);
for i = 1 : N_delta
    Parameters.delta = delta_vector(i);
    disp(['Progress: (i,j) = (', num2str(i), ',:) out of (', num2str(N_delta), ',', num2str(N_theta), ')']);
    parfor j = 1 : N_theta
        localParameters = Parameters;
        localParameters.theta = theta_vector(j);
        [StationaryMeanX, StationaryVarianceX, ~] = ComputeStatistics_Stationary(PropensityFunction, StoichiometryMatrix, localParameters, IC, TimeSpan_Long, TransientPortion, N_Trajectories_Few);
        StationaryMean(i,j) = StationaryMeanX(Output_Index);
        StationaryVariance(i,j) = StationaryVarianceX(Output_Index);
    end
end
SimTime = toc(SimTime);
disp(['Successfully completed the run of: "', mfilename, '.m". Simulation Time is: ', num2str(SimTime/3600), ' hours.']);

%% Save data
Save_Flag = 1;
FileName = fullfile(DataPath, 'FP_ThreeSpecies_delta_theta');
if Save_Flag == 1
    clear euler pool
    save(FileName);
end