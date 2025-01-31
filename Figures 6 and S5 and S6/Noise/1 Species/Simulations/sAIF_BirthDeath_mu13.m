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
StoichiometryMatrix = StoichiometryMatrix_sAIF_BirthDeath();
PropensityFunction = @PropensityFunction_sAIF_BirthDeath;
Parameters = Parameters_sAIF_BirthDeath();
Output_Index = 1;
N_mu = 32;
mu_vector = linspace(0.1, 10, N_mu);
Script_Index = 13;
Parameters.mu = mu_vector(Script_Index);

%% Simulation Settings
IC = [0; 0; 0];
TimeSpan_Long = [0, 1e6];
N_Trajectories_Few = 1;
TransientPortion = 0.01;

%% Swept Parameters
N_delta = 16;
N_theta = 32;
N_eta = 32;
delta_vector = logspace(-3, -1, N_delta);
theta_vector = logspace(-5, 1, N_theta);
eta_vector = logspace(-5, 5, N_eta);

%% Simulations
SimTime = tic;
StationaryMean = zeros(N_delta, N_theta, N_eta);
StationaryVariance = zeros(N_delta, N_theta, N_eta);
for i = 1 : N_delta
    Parameters.delta = delta_vector(i);
    for j = 1 : N_theta
        Parameters.theta = theta_vector(j);
        disp(['Progress: (i,j,k) = (', num2str(i), ',' num2str(j), ',:) out of (', num2str(N_delta), ',', num2str(N_theta), ',', num2str(N_eta), ')']);
        parfor k = 1 : N_eta
            localParameters = Parameters;
            localParameters.eta = eta_vector(k);
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
FileName = fullfile(DataPath, ['sAIF_BirthDeath_delta_theta_eta_mu', num2str(Script_Index)]);
if Save_Flag == 1
    clear euler pool
    save(FileName);
end