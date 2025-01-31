%% Clear Workspace
close all
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
FixedPoint = @ComputeFP_sAIF_BirthDeath;
CheckStability = @ComputeStability_sAIF_BirthDeath;
Output_Index = 1;

%% Fixed Parameters
Parameters.gamma_1 = 0.1;
Parameters.delta = 0.1;
Parameters.alpha = 2;
Parameters.kappa = 0.05;

%% Disturbance
Disturbance = 'gamma_1';
Disturbance_Factor = 1/2;

%% Grid
N_mu = 300;
N_eta = 60;
N_theta = 300;
mu_vector = linspace(0, 10, N_mu);
eta_vector = [0, logspace(-3, 5, N_eta-1)];
theta_vector = logspace(-6, 1, N_theta);

%% Simulations for sAIF
SimTime = tic;
SS_PreDisturbance = zeros(N_mu, N_eta, N_theta);
SS_PostDisturbance = zeros(N_mu, N_eta, N_theta);
Stability_PreDisturbance = zeros(N_mu, N_eta, N_theta);
Stability_PostDisturbance = zeros(N_mu, N_eta, N_theta);
for i = 1 : N_mu
    Parameters.mu = mu_vector(i);
    for j = 1 : N_eta
        Parameters.eta = eta_vector(j);
        disp(['Progress: (i,j,k) = (', num2str(i), ',', num2str(j), ',:) out of (', num2str(N_mu), ',', num2str(N_eta), ',', num2str(N_theta), ')']);
        for k = 1 : N_theta
            localParameters = Parameters;
            localParameters.theta = theta_vector(k);
            DisturbedParameters = localParameters;
            DisturbedParameters.(Disturbance) = localParameters.(Disturbance) * Disturbance_Factor;
            FP = FixedPoint(localParameters); 
            SS_PreDisturbance(i,j,k) = FP(Output_Index);
            Stability_PreDisturbance(i,j,k) = CheckStability(localParameters, FP);
            FP = FixedPoint(DisturbedParameters); 
            SS_PostDisturbance(i,j,k) = FP(Output_Index);
            Stability_PostDisturbance(i,j,k) = CheckStability(DisturbedParameters, FP);
        end
    end
end
SimTime = toc(SimTime);
SS_Error = abs(SS_PreDisturbance - SS_PostDisturbance) ./ SS_PreDisturbance;
disp(['Successfully completed the run of: "', mfilename, '.m". Simulation Time is: ', num2str(SimTime/3600), ' hours.']);

%% Infinite eta
SS_PreDisturbance_Inf = zeros(N_mu, N_theta);
SS_PostDisturbance_Inf = zeros(N_mu, N_theta);
for i = 1 : N_mu
    Parameters.mu = mu_vector(i);
    for j = 1 : N_theta
        localParameters = Parameters;
        localParameters.theta = theta_vector(j);
        DisturbedParameters = localParameters;
        DisturbedParameters.(Disturbance) = localParameters.(Disturbance) * Disturbance_Factor;
        SS_PreDisturbance_Inf(i,j) = ComputeFP_sAIF_BirthDeath_Inf_eta(localParameters); 
        SS_PostDisturbance_Inf(i,j) = ComputeFP_sAIF_BirthDeath_Inf_eta(DisturbedParameters); 
    end
end
SS_Error_Inf = abs(SS_PreDisturbance_Inf - SS_PostDisturbance_Inf) ./ SS_PreDisturbance_Inf;

%% Contour for Fixed theta
Setpoint_Fixed_theta = 5;
N_Contour = 100;
eta_Contour = [1e-10, logspace(-3, 5, N_Contour-1)];
z_2 = Parameters.kappa * (Parameters.alpha/Parameters.gamma_1/Setpoint_Fixed_theta - 1);
theta_Fixed = (Parameters.delta * Setpoint_Fixed_theta) * z_2;
mu_Contour = ((eta_Contour * z_2 + Parameters.delta) ./ eta_Contour) * (theta_Fixed * Setpoint_Fixed_theta / z_2 - Parameters.delta);
DisturbedParameters_Fixed_theta = Parameters;
DisturbedParameters_Fixed_theta.theta = theta_Fixed;
DisturbedParameters_Fixed_theta.(Disturbance) = Parameters.(Disturbance) * Disturbance_Factor;
SS_Error_Fixed_theta = zeros(N_Contour, 1);
for i = 1 : length(eta_Contour)
    DisturbedParameters_Fixed_theta.mu = mu_Contour(i);
    DisturbedParameters_Fixed_theta.eta = eta_Contour(i);
    FP = FixedPoint(DisturbedParameters_Fixed_theta);
    SS_PostDisturbance_Fixed_theta = FP(Output_Index);
    SS_Error_Fixed_theta(i) = abs(Setpoint_Fixed_theta - SS_PostDisturbance_Fixed_theta) ./ Setpoint_Fixed_theta;
end

%% Contour for Fixed mu
Setpoint_Fixed_mu = 10;
mu_Fixed = 5;
z_2 = Parameters.kappa * (Parameters.alpha/Parameters.gamma_1/Setpoint_Fixed_mu - 1);
theta_Contour = (mu_Fixed * eta_Contour ./ (eta_Contour * z_2 + Parameters.delta) + Parameters.delta) * (z_2/Setpoint_Fixed_mu);
DisturbedParameters_Fixed_mu = Parameters;
DisturbedParameters_Fixed_mu.mu = mu_Fixed;
DisturbedParameters_Fixed_mu.(Disturbance) = Parameters.(Disturbance) * Disturbance_Factor;
SS_Error_Fixed_mu = zeros(N_Contour, 1);
for i = 1 : length(eta_Contour)
    DisturbedParameters_Fixed_mu.theta = theta_Contour(i);
    DisturbedParameters_Fixed_mu.eta = eta_Contour(i);
    FP = FixedPoint(DisturbedParameters_Fixed_mu);
    SS_PostDisturbance_Fixed_mu = FP(Output_Index);
    SS_Error_Fixed_mu(i) = abs(Setpoint_Fixed_mu - SS_PostDisturbance_Fixed_mu) ./ Setpoint_Fixed_mu;
end

%% Save data
Save_Flag = 1;
FileName = fullfile(DataPath, 'sAIF_BirthDeath_SSError_mu_eta_theta');
if Save_Flag == 1
    clear euler pool
    save(FileName);
end