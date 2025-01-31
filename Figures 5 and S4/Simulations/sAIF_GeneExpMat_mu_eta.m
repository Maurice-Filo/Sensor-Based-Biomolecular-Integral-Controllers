%% Clear Workspace
close all
clear
clc
disp(['Running: "', mfilename, '.m"']);

%% Configure Working Environment
ProjectPath = '/cluster/home/mfilo/Noise Ideal';
DataPath = fullfile(ProjectPath, 'Data');
addpath(genpath(ProjectPath));
euler = parcluster('local');
pool = parpool(euler, 16);

% ProjectPath = '/Users/mfilo/Library/CloudStorage/GoogleDrive-maurice.g.filo@gmail.com/My Drive/Research at ETH/My Papers and Presentations/Bacterial PI/Matlab/Noise/Noise Ideal';
% DataPath = fullfile(ProjectPath, 'Data');
% addpath(genpath(ProjectPath));

%% Network Construction
StoichiometryMatrix = StoichiometryMatrix_sAIF_GeneExpMature();
PropensityFunction = @PropensityFunction_sAIF_GeneExpMature;
Parameters = Parameters_sAIF_GeneExpMature();
Output_Index = 3;

%% Simulation Settings
IC = [0; 0; 0; 0; 0];
TimeSpan_Long = [0, 1e6];
N_Trajectories_Few = 1;
TransientPortion = 0.01;

%% Swept Parameters
N_mu = 48;
N_eta = 24;
mu_vector = linspace(1, 10, N_mu);
eta_vector = logspace(-2, 2, N_eta);

%% Simulations
SimTime = tic;
StationaryMean = zeros(N_mu, N_eta);
StationaryVariance = zeros(N_mu, N_eta);
for j = 1 : N_eta
    Parameters.eta = eta_vector(j);
    disp(['Progress: (i,j) = (', num2str(j), ',:) out of (', num2str(N_eta), ',' num2str(N_mu), ')']);
    parfor i = 1 : N_mu
        localParameters = Parameters;
        localParameters.mu = mu_vector(i);
        [StationaryMeanX, StationaryVarianceX, ~] = ComputeStatistics_Stationary(PropensityFunction, StoichiometryMatrix, localParameters, IC, TimeSpan_Long, TransientPortion, N_Trajectories_Few);
        StationaryMean(i,j) = StationaryMeanX(Output_Index);
        StationaryVariance(i,j) = StationaryVarianceX(Output_Index);
    end
end
SimTime = toc(SimTime);
disp(['Successfully completed the run of: "', mfilename, '.m". Simulation Time is: ', num2str(SimTime/3600), ' hours.']);

%% Save data
Save_Flag = 1;
FileName = fullfile(DataPath, 'sAIF_GeneExpMat_mu_eta');
if Save_Flag == 1
    clear euler pool
    save(FileName);
end