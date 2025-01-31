%% Clear Workspace
clear;
close all;
clc;
Save_Flag = 0;

%% Load Data
Data_sAIF = load('sAIF_BirthDeath_SSError_mu_eta_theta.mat');
SS_PreDisturbance = Data_sAIF.SS_PreDisturbance;
Stability_PreDisturbance = Data_sAIF.Stability_PreDisturbance;
Stability_PostDisturbance = Data_sAIF.Stability_PostDisturbance;
SS_Error = Data_sAIF.SS_Error;
N_mu = Data_sAIF.N_mu;
N_eta = Data_sAIF.N_eta;
N_theta = Data_sAIF.N_theta;
mu_vector = Data_sAIF.mu_vector;
eta_vector = Data_sAIF.eta_vector;
theta_vector = Data_sAIF.theta_vector;
Output_Index = Data_sAIF.Output_Index;
Parameters = Data_sAIF.Parameters;
Disturbance = Data_sAIF.Disturbance;
Disturbance_Factor = Data_sAIF.Disturbance_Factor;

% Infinite eta
SS_PreDisturbance_Inf = Data_sAIF.SS_PreDisturbance_Inf;
Error_Inf = Data_sAIF.SS_Error_Inf;

% Contour for Fixed mu
SS_Error_Fixed_mu = Data_sAIF.SS_Error_Fixed_mu;
Setpoint_Fixed_mu = Data_sAIF.Setpoint_Fixed_mu;
eta_Contour = Data_sAIF.eta_Contour;
theta_Contour = Data_sAIF.theta_Contour;
mu_Fixed = Data_sAIF.mu_Fixed;

%% Selected Simulations for Fixed mu
selected_indices_Fixed_mu = [1, 51, 100];
eta_selected_Fixed_mu = eta_Contour(selected_indices_Fixed_mu); 
theta_selected_Fixed_mu = theta_Contour(selected_indices_Fixed_mu);
SS_Error_Fixed_mu_selected = SS_Error_Fixed_mu(selected_indices_Fixed_mu);

% Simulation Settings
IC = zeros(3,1);
tf = 200;
Nt = 1000;
Solver = 'ODE15s';

% Simulations
y_Fixed_mu = zeros(length(eta_selected_Fixed_mu), 2*Nt-1);
Parameters.mu = mu_Fixed;
for i = 1 : length(eta_selected_Fixed_mu)
    Parameters.theta = theta_selected_Fixed_mu(i);
    Parameters.eta = eta_selected_Fixed_mu(i);
    DisturbedParameters = Parameters;
    DisturbedParameters.(Disturbance) = Parameters.(Disturbance) * Disturbance_Factor;
    [time_vector1, x1] = DSA(Data_sAIF.StoichiometryMatrix, Data_sAIF.PropensityFunction, Parameters, IC, tf, Nt, Solver);
    [time_vector2, x2] = DSA(Data_sAIF.StoichiometryMatrix, Data_sAIF.PropensityFunction, DisturbedParameters, x1(:,end), tf, Nt, Solver);
    time_vector = [time_vector1, time_vector2(2:end) + time_vector1(end)];
    x = [x1, x2(:,2:end)];
    y_Fixed_mu(i,:) = x(Output_Index,:);
end

%% Figure Settings
SS = 4;
Figure_Width = 7.5 * SS;
Figure_Height = 2 * SS;
FontSize = 5 * SS;
FontSize_Small = 3 * SS;
FontSize_Large = 6 * SS;
LineWidth = 4;
LineWidth_Thin = 0.00001;
MarkerSize = 10;
NominalColors = [55,126,184; ... % Blue
                 77,175,74; ... % Green
                 255,127,0; ... % Orange
                 228,26,28; ... % Red
                 163, 124, 182; ... % Purple
                 255, 0, 255]/255; ... % Magenta

%% Set Figure 1
Figure1_Name = 'SS_Error_BirthDeath_mu_eta_theta2D';
Handle_Figure1 = figure;
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [0, Figure_Height, Figure_Width/2.2, Figure_Height];
    Handle_Figure1.PaperPositionMode = 'auto';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];

%% Set Axis 1: Steady-State Error
Handle_Axis1 = axes(Handle_Figure1);
    Handle_Axis1.Position = [0.11, 0.13, 0.85, 0.82];
    Handle_Axis1.Box = 'on';
    Handle_Axis1.FontSize = FontSize;
    hold(Handle_Axis1, 'on');
    grid(Handle_Axis1, 'on');
    Handle_Axis1.Layer = 'top';
    Handle_Axis1.XMinorGrid = 'off';
    Handle_Axis1.YMinorGrid = 'off';
    Handle_Axis1.XScale = 'linear';
    Handle_Axis1.YScale = 'linear';
    Handle_Axis1.XLabel.Interpreter = 'latex';
    Handle_Axis1.YLabel.Interpreter = 'latex';
    Handle_Axis1.XLim = [0, 20];
    Handle_Axis1.YLim = [0.95*min(min(min(Error_Inf))), 1];
    Handle_Axis1.YTick = [1e-3, 1e-2, 1e-1, 1];
    % Handle_Axis1.XLabel.String = 'Setpoint';
    % Handle_Axis1.YLabel.String = 'Error';
    Handle_Axis1.YScale = 'Log';

% Scatter Points
Indices_mu = 1 :1: N_mu;
Indices_eta = 1 :1: N_eta;
Indices_theta = 1 :1: N_theta;
SetPoints = SS_PreDisturbance(Indices_mu, Indices_eta, Indices_theta); SetPoints = SetPoints(:);
Errors = SS_Error(Indices_mu, Indices_eta, Indices_theta);  Errors = Errors(:);
scatter(Handle_Axis1, SetPoints, Errors, MarkerSize, 'MarkerFaceColor', 'flat', 'CData', NominalColors(1,:), 'Marker', 'o', 'MarkerEdgeAlpha', 0.1, 'MarkerFaceAlpha', 1, 'MarkerEdgeColor', 'k', 'LineWidth', LineWidth_Thin);

% N_x = 60; N_y = 40;
% x_Grid = linspace(0, 20, N_x);
% y_Grid = logspace(-3, 0, N_y);
% [x_Grid, y_Grid] = meshgrid(x_Grid, y_Grid);
% Selected_x = zeros(numel(x_Grid),1);
% Selected_y = zeros(numel(x_Grid),1);
% for i = 1 : numel(x_Grid)
%     Grid_Point = [x_Grid(i), y_Grid(i)];
%     Distances = sqrt( (Grid_Point(1) - SetPoints).^2 + (Grid_Point(2) - Errors).^2);
%     [~, Index] = min(Distances);
%     Selected_x(i) = SetPoints(Index);
%     Selected_y(i) = Errors(Index);
% end
% scatter(Handle_Axis1, Selected_x, Selected_y, MarkerSize, 'MarkerFaceColor', 'flat', 'CData', NominalColors(1,:), 'Marker', 'o', 'MarkerEdgeAlpha', 1, 'MarkerFaceAlpha', 1, 'MarkerEdgeColor', 'k', 'LineWidth', LineWidth_Thin);

% Sandwiching Curves
SetPoints = squeeze(SS_PreDisturbance(1,1,:));  
Errors = squeeze(SS_Error(1,1,:)); 
plot(Handle_Axis1, SetPoints, Errors, 'Color', NominalColors(5,:), 'LineWidth', LineWidth*1.5);

% Simulation Points
Colors = [NominalColors(5,:); NominalColors(3,:); NominalColors(2,:)];
scatter(Handle_Axis1, 10*ones(length(SS_Error_Fixed_mu_selected), 1), SS_Error_Fixed_mu_selected, MarkerSize*10, 'MarkerFaceColor', 'flat', 'LineWidth', LineWidth_Thin, 'CData', Colors, 'MarkerEdgeColor', 'k', 'LineWidth', LineWidth_Thin);

%% Set Figure 2
Figure2_Name = 'SS_Error_BirthDeath_Simulations';
Handle_Figure2 = figure;
    Handle_Figure2.PaperUnits = 'centimeters';
    Handle_Figure2.Units = 'centimeters';
    Handle_Figure2.Position = [0, Figure_Height, Figure_Width/2, Figure_Height/1.5];
    Handle_Figure2.PaperPositionMode = 'auto';
    Handle_Figure2.PaperSize = [Handle_Figure2.PaperPosition(3), Handle_Figure2.PaperPosition(4)];

%% Set Axis 2
Handle_Axis2 = axes(Handle_Figure2);
    Handle_Axis2.Position = [0.07, 0.19, 0.89, 0.75];
    Handle_Axis2.Box = 'on';
    Handle_Axis2.FontSize = FontSize;
    hold(Handle_Axis2, 'on');
    grid(Handle_Axis2, 'on');
    Handle_Axis2.Layer = 'top';
    Handle_Axis2.XMinorGrid = 'off';
    Handle_Axis2.YMinorGrid = 'off';
    Handle_Axis2.XScale = 'linear';
    Handle_Axis2.YScale = 'linear';
    Handle_Axis2.XLabel.Interpreter = 'latex';
    Handle_Axis2.YLabel.Interpreter = 'latex';
    Handle_Axis2.XLim = [150, 300];
    Handle_Axis2.XTick = 150:50:300;
    Handle_Axis2.XTickLabel = {'', '0', '50', '100'};
    Handle_Axis2.YLim = [5, 20];
    Handle_Axis2.YTick = [5, 10, 15, 20];
Colors = [NominalColors(5,:); NominalColors(3,:); NominalColors(2,:)];
for i = 1 : length(eta_selected_Fixed_mu)
    plot(Handle_Axis2, time_vector, y_Fixed_mu(i,:), 'LineWidth', LineWidth / 2, 'Color', Colors(i,:));
end
% Handle_Legend = legend(Handle_Axis2, {'filtered P (eta = 0)', ['non-ideal sAIF (\eta = ', num2str(eta_selected_Fixed_mu(2)), ')'], 'non-ideal sAIF (\eta = 10^5)'});
% Handle_Legend.Location = 'best';
% Handle_Legend.FontSize = FontSize/2;

%% Save Figures
if Save_Flag == 1  
    Handle_Figure1.Color = [1, 1, 1];
    % Handle_Figure1.Color = 'none';
    % set(Handle_Figure1, 'InvertHardCopy', 'off');
    print(Handle_Figure1, Figure1_Name, '-dpng'); 
    % Handle_Figure1.Color = [1, 1, 1];
    
    Handle_Figure2.Color = [1, 1, 1];
    Handle_Figure2.Color = 'none';
    set(Handle_Figure2, 'InvertHardCopy', 'off');
    print(Handle_Figure2, Figure2_Name, '-dpdf', '-vector','-bestfit');
    Handle_Figure2.Color = [1, 1, 1];
end

