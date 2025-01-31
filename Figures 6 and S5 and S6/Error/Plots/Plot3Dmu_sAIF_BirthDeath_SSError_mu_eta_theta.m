%% Clear Workspace
clear;
close all;
clc;
Save_Flag = 0;

%% Load Data
% Finite eta
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

% Contour for Fixed theta
SS_Error_Fixed_theta = Data_sAIF.SS_Error_Fixed_theta;
Setpoint_Fixed_theta = Data_sAIF.Setpoint_Fixed_theta;
eta_Contour = Data_sAIF.eta_Contour;
mu_Contour = Data_sAIF.mu_Contour;
theta_Fixed = Data_sAIF.theta_Fixed;

% Contour for Fixed mu
SS_Error_Fixed_mu = Data_sAIF.SS_Error_Fixed_mu;
Setpoint_Fixed_mu = Data_sAIF.Setpoint_Fixed_mu;
theta_Contour = Data_sAIF.theta_Contour;
mu_Fixed = Data_sAIF.mu_Fixed;

%% Selected Simulations for Fixed theta
selected_indices_Fixed_theta = [1, 20, 25, 30, 35, 40, 45, 100];
eta_selected_Fixed_theta = eta_Contour(selected_indices_Fixed_theta); 
mu_selected_Fixed_theta = mu_Contour(selected_indices_Fixed_theta);
SS_Error_Fixed_theta_selected = SS_Error_Fixed_theta(selected_indices_Fixed_theta);

% Simulation Settings
IC = zeros(3,1);
tf = 100;
Nt = 1000;
Solver = 'ODE15s';

% Simulations
y_Fixed_theta = zeros(length(eta_selected_Fixed_theta), 2*Nt-1);
Parameters.theta = theta_Fixed;
for i = 1 : length(eta_selected_Fixed_theta)
    Parameters.mu = mu_selected_Fixed_theta(i);
    Parameters.eta = eta_selected_Fixed_theta(i);
    DisturbedParameters = Parameters;
    DisturbedParameters.(Disturbance) = Parameters.(Disturbance) * Disturbance_Factor;
    [time_vector1, x1] = DSA(Data_sAIF.StoichiometryMatrix, Data_sAIF.PropensityFunction, Parameters, IC, tf, Nt, Solver);
    [time_vector2, x2] = DSA(Data_sAIF.StoichiometryMatrix, Data_sAIF.PropensityFunction, DisturbedParameters, x1(:,end), tf, Nt, Solver);
    time_vector = [time_vector1, time_vector2(2:end) + time_vector1(end)];
    x = [x1, x2(:,2:end)];
    y_Fixed_theta(i,:) = x(Output_Index,:);
end

%% Selected Simulations for Fixed mu
selected_indices_Fixed_mu = [1, 35, 40, 45, 50, 55, 60, 100];
eta_selected_Fixed_mu = eta_Contour(selected_indices_Fixed_mu); 
theta_selected_Fixed_mu = theta_Contour(selected_indices_Fixed_mu);
SS_Error_Fixed_mu_selected = SS_Error_Fixed_mu(selected_indices_Fixed_mu);

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
Figure_Height = 3 * SS;
FontSize = 5 * SS;
FontSize_Small = 3 * SS;
FontSize_Large = 6 * SS;
LineWidth = 4;
LineWidth_Thin = 0.01;
MarkerSize = 5;
NominalColors = [55,126,184; ... % Blue
                 77,175,74; ... % Green
                 255,0,255; ... % Magenta
                 228,26,28; ... % Red
                 127, 0, 255]/255; ... % Purple

%% Set Figure 1
Figure1_Name = 'SS_Error_BirthDeath_mu_eta_theta';
Handle_Figure1 = figure;
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [0, Figure_Height, Figure_Width/2, Figure_Height];
    Handle_Figure1.PaperPositionMode = 'auto';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];

%% Set Axis 1: Steady-State Error
Handle_Axis1 = axes(Handle_Figure1);
    Handle_Axis1.Position = [0.12, 0.1, 0.65, 0.85];
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
    Handle_Axis1.XLim = [0, mu_vector(end)];
    Handle_Axis1.YLim = [0, 20];
    Handle_Axis1.ZLim = [0.95*min(min(min(Error_Inf))), 1];
    Handle_Axis1.ZTick = [1e-3, 1e-2, 1e-1, 1];
    % Handle_Axis1.XLabel.String = 'mu';
    % Handle_Axis1.YLabel.String = 'Setpoint';
    % Handle_Axis1.ZLabel.String = 'Error';
    Handle_Axis1.ZScale = 'Log';

% Scatter Points
Indices_mu = 1 :2: N_mu;
Indices_eta = 1 :1: N_eta;
Indices_theta = 1 :2: N_theta;
cmap = colormap(turbo(length(Indices_eta)));
% ColorIndices = reshape(repmat(1:length(Indices_eta), length(Indices_theta), 1, length(Indices_mu)), [], 1); 
[~, ColorIndices, ~] = ndgrid(1:length(Indices_theta), 1:length(Indices_eta), 1:length(Indices_mu));
ColorIndices = ColorIndices(:);
Colors = cmap(ColorIndices, :);
SetPoints = SS_PreDisturbance(Indices_mu, Indices_eta, Indices_theta); SetPoints = SetPoints(:);
Errors = SS_Error(Indices_mu, Indices_eta, Indices_theta);  Errors = Errors(:);
[Mus, ~, ~] = ndgrid(mu_vector(Indices_mu), eta_vector(Indices_eta), theta_vector(Indices_theta));
Mus = Mus(:);
scatter3(Handle_Axis1, Mus, SetPoints, Errors, MarkerSize, 'MarkerFaceColor', 'flat', 'CData', Colors, 'Marker', 'o', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 0);
clim([eta_vector(1), eta_vector(end)]);
Handle_Colorbar = colorbar();
Handle_Colorbar.Location = 'manual';
Handle_Colorbar.Position = [0.8553, 0.03, 0.0376, 0.88];
Handle_Axis1.View = [31.5241, 16.4699];
Handle_Colorbar.Ruler.TickValues = [min(eta_vector(Indices_eta)), max(eta_vector(Indices_eta))]; 

% Sandwiching Surfaces
[Mus, ~] = meshgrid(mu_vector, theta_vector);
SetPoints = squeeze(SS_PreDisturbance(:,1,:));  
Errors = squeeze(SS_Error(:,1,:)); 
Handle_Surface = surf(Handle_Axis1, Mus', SetPoints, Errors);
Handle_Surface.FaceColor = cmap(1,:);
Handle_Surface.FaceAlpha = 0.8;
Handle_Surface.EdgeColor = 'none';

SetPoints = SS_PreDisturbance_Inf;
Errors = Error_Inf;
[xq, yq] = meshgrid(linspace(min(Mus(:)), max(Mus(:)), 600), linspace(min(SetPoints(:)), max(SetPoints(:)), 600));
zq = griddata(Mus', SetPoints, Errors, xq, yq, 'linear');  
Handle_Surface = surf(Handle_Axis1, xq, yq, zq);
Handle_Surface.FaceColor = cmap(end,:);
Handle_Surface.FaceAlpha = 0.8;
Handle_Surface.EdgeColor = 'none';

%% Set Figure 2
Figure2_Name = 'SS_Error_BirthDeath_mu_eta_theta_Simulations';
Handle_Figure2 = figure;
    Handle_Figure2.PaperUnits = 'centimeters';
    Handle_Figure2.Units = 'centimeters';
    Handle_Figure2.Position = [0, Figure_Height, Figure_Width/2, Figure_Height];
    Handle_Figure2.PaperPositionMode = 'auto';
    Handle_Figure2.PaperSize = [Handle_Figure2.PaperPosition(3), Handle_Figure2.PaperPosition(4)];

%% Set Axis 2
Handle_Axis2 = axes(Handle_Figure2);
    Handle_Axis2.Position = [0.07, 0.1, 0.89, 0.75];
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
for i = 1 : length(eta_selected_Fixed_theta)
    plot(Handle_Axis2, time_vector, y_Fixed_theta(i,:), 'LineWidth', LineWidth / 2, 'Color', [NominalColors(3,:), 0.1 + (i-1) * (0.7 - 0.1)/(length(eta_selected_Fixed_theta) - 1)]);
end

for i = 1 : length(eta_selected_Fixed_mu)
    plot(Handle_Axis2, time_vector, y_Fixed_mu(i,:), 'LineWidth', LineWidth / 2, 'Color', [NominalColors(1,:), 0.1 + (i-1) * (0.7 - 0.1)/(length(eta_selected_Fixed_mu) - 1)]);
end

Handle_Axis3 = axes(Handle_Figure2);
    Handle_Axis3.Position = [0.18, 0.6, 0.35, 0.4];
    Handle_Axis3.Box = Handle_Axis1.Box;
    hold(Handle_Axis3, 'on');
    grid(Handle_Axis3, 'on');
    Handle_Axis3.Layer = Handle_Axis1.Layer;
    Handle_Axis3.XMinorGrid = Handle_Axis1.XMinorGrid;
    Handle_Axis3.YMinorGrid = Handle_Axis1.YMinorGrid;
    Handle_Axis3.XScale = Handle_Axis1.XScale;
    Handle_Axis3.YScale = Handle_Axis1.YScale;
    Handle_Axis3.XLim = Handle_Axis1.XLim;
    Handle_Axis3.YLim = Handle_Axis1.YLim;
    Handle_Axis3.ZLim = Handle_Axis1.ZLim;
    Handle_Axis3.XTick = [];
    Handle_Axis3.YTick = [];
    Handle_Axis3.ZTick = [];
    Handle_Axis3.ZScale = Handle_Axis1.ZScale;
    Handle_Axis3.View = Handle_Axis1.View;

% Scatter Points
SetPoints = reshape(SS_PreDisturbance(Indices_mu, Indices_eta, Indices_theta), [], 1);  
Errors = reshape(SS_Error(Indices_mu, Indices_eta, Indices_theta), [], 1);  
Mus = repmat(mu_vector(Indices_mu)', length(Indices_theta) * length(Indices_eta), 1);
scatter3(Handle_Axis3, Mus, SetPoints, Errors, MarkerSize, 'MarkerFaceColor', 'flat', 'CData', Colors, 'Marker', 'o', 'MarkerFaceAlpha', 0.03, 'MarkerEdgeAlpha', 0);

% Sandwiching Surfaces
[Mus, ~] = meshgrid(mu_vector, theta_vector);
SetPoints = squeeze(SS_PreDisturbance(:,1,:));  
Errors = squeeze(SS_Error(:,1,:)); 
Handle_Surface = surf(Handle_Axis3, Mus', SetPoints, Errors);
Handle_Surface.FaceColor = cmap(1,:);
Handle_Surface.FaceAlpha = 0.2;
Handle_Surface.EdgeColor = 'none';

SetPoints = SS_PreDisturbance_Inf;
Errors = Error_Inf;
[xq, yq] = meshgrid(linspace(min(Mus(:)), max(Mus(:)), 1000), linspace(min(SetPoints(:)), max(SetPoints(:)), 1000));
zq = griddata(Mus', SetPoints, Errors, xq, yq, 'linear');  
Handle_Surface = surf(Handle_Axis3, xq, yq, zq);
Handle_Surface.FaceColor = cmap(end,:);
Handle_Surface.FaceAlpha = 0.2;
Handle_Surface.EdgeColor = 'none';

% Lines for Fixed theta
plot3(Handle_Axis3, mu_Contour, Setpoint_Fixed_theta * ones(length(mu_Contour), 1), SS_Error_Fixed_theta, 'Color', NominalColors(3,:), 'LineWidth', LineWidth);
% plot3(Handle_Axis3, mu_selected_Fixed_theta, Setpoint_Fixed_theta*ones(length(mu_selected_Fixed_theta), 1), SS_Error_Fixed_theta_selected, 'Marker', 's', 'LineStyle', 'none', 'MarkerSize', MarkerSize*3, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
PlotCone(Handle_Axis3, mu_Contour, Setpoint_Fixed_theta * ones(length(mu_Contour), 1), SS_Error_Fixed_theta, 27, 0, 2, 0.2, 2, NominalColors(3,:));

% Lines for Fixed mu
plot3(Handle_Axis3, mu_Fixed * ones(length(theta_Contour), 1), Setpoint_Fixed_mu * ones(length(theta_Contour), 1), SS_Error_Fixed_mu, 'Color', NominalColors(1,:), 'LineWidth', LineWidth);
% plot3(Handle_Axis3, mu_Fixed * ones(length(theta_selected_Fixed_mu), 1), Setpoint_Fixed_mu*ones(length(theta_selected_Fixed_mu), 1), SS_Error_Fixed_mu_selected, 'Marker', 's', 'LineStyle', 'none', 'MarkerSize', MarkerSize*3, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
PlotCone(Handle_Axis3, mu_Fixed * ones(length(theta_Contour), 1), Setpoint_Fixed_mu * ones(length(theta_Contour), 1), SS_Error_Fixed_mu, 70, 0, 0.02, 1, 1, NominalColors(1,:));

%% Save Figures
if Save_Flag == 1  
    Handle_Figure1.Color = [1, 1, 1];
    % Handle_Figure1.Color = 'none';
    set(Handle_Figure1, 'InvertHardCopy', 'off');
    print(Handle_Figure1, Figure1_Name, '-dpng', '-r600'); 
    
    Handle_Figure2.Color = [1, 1, 1];
    % Handle_Figure2.Color = 'none';
    set(Handle_Figure2, 'InvertHardCopy', 'off');
    print(Handle_Figure2, Figure2_Name, '-dpng', '-r600');
end

