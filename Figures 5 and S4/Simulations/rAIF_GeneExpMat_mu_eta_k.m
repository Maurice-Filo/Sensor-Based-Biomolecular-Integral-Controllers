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
StoichiometryMatrix = StoichiometryMatrix_rAIF_GeneExpMature();
PropensityFunction = @PropensityFunction_rAIF_GeneExpMature;
Parameters = Parameters_rAIF_GeneExpMature();
Output_Index = 3;

%% Simulation Settings
IC = [0; 0; 0; 0; 0];
TimeSpan_Long = [0, 1e6];
N_Trajectories_Few = 1;
TransientPortion = 0.01;

%% Swept Parameters
N_mu = 48;
N_eta = 24;
N_k = 12;
mu_vector = linspace(1, 10, N_mu);
eta_vector = logspace(-2, 2, N_eta);
k_vector = linspace(1e-3, 1, N_k);

%% Simulations
SimTime = tic;
StationaryMean = zeros(N_mu, N_eta, N_k);
StationaryVariance = zeros(N_mu, N_eta, N_k);
for k = 1 : N_k
    Parameters.k = k_vector(k);
    for j = 1 : N_eta
        Parameters.eta = eta_vector(j);
        disp(['Progress: (i,j,k) = (', num2str(k), ',', num2str(j), ',:) out of (', num2str(N_k), ',' num2str(N_eta), ',', num2str(N_mu), ')']);
        parfor i = 1 : N_mu
            localParameters = Parameters;
            localParameters.mu = mu_vector(i);
            [StationaryMeanX, StationaryVarianceX, ~] = ComputeStatistics_Stationary(PropensityFunction, StoichiometryMatrix, localParameters, IC, TimeSpan_Long, TransientPortion, N_Trajectories_Few);
            StationaryMean(i,j,k) = StationaryMeanX(Output_Index);
            StationaryVariance(i,j,k) = StationaryVarianceX(Output_Index);
        end
    end
end
SimTime = toc(SimTime);
disp(['Successfully completed the run of: "', mfilename, '.m". Simulation Time is: ', num2str(SimTime/3600), ' hours.']);

%% Save data
Save_Flag = 1;
FileName = fullfile(DataPath, 'rAIF_GeneExpMat_mu_eta_k');
if Save_Flag == 1
    clear euler pool
    save(FileName);
end