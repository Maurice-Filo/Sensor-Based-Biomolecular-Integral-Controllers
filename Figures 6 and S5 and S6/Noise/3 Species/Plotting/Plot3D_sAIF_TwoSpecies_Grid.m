%% Clear Workspace
clear;
close all;
clc;
Save_Flag = 0;
surface = true;

%% Load and Process sAIF Data
Data_sAIF_mu_1 = load('sAIF_TwoSpecies_delta_theta_eta_mu1.mat');
delta_vector_sAIF = Data_sAIF_mu_1.delta_vector;
theta_vector_sAIF = Data_sAIF_mu_1.theta_vector;
eta_vector_sAIF = Data_sAIF_mu_1.eta_vector;
mu_vector_sAIF = Data_sAIF_mu_1.mu_vector;
N_delta_sAIF = Data_sAIF_mu_1.N_delta;
N_theta_sAIF = Data_sAIF_mu_1.N_theta;
N_eta_sAIF = Data_sAIF_mu_1.N_eta;
N_mu_sAIF = Data_sAIF_mu_1.N_mu;
Means_sAIF = zeros(N_delta_sAIF, N_theta_sAIF, N_eta_sAIF, N_mu_sAIF);
CVs_sAIF = zeros(N_delta_sAIF, N_theta_sAIF, N_eta_sAIF, N_mu_sAIF);
for i = 1 : N_mu_sAIF
    Data_sAIF_mu_i = load(['sAIF_TwoSpecies_delta_theta_eta_mu', num2str(i), '.mat']);
    Means_sAIF(:,:,:,i) = Data_sAIF_mu_i.StationaryMean;
    CVs_sAIF(:,:,:,i) = sqrt(Data_sAIF_mu_i.StationaryVariance) ./ Data_sAIF_mu_i.StationaryMean;
end

% Filter Data
delta_Selectors_sAIF = 1:N_delta_sAIF;
theta_Selectors_sAIF = 1:N_theta_sAIF;
eta_Selectors_sAIF = 1:N_eta_sAIF;
mu_Selectors_sAIF = 1:N_mu_sAIF;

delta_Masks_sAIF = ismember(1:N_delta_sAIF, delta_Selectors_sAIF);
theta_Masks_sAIF = ismember(1:N_theta_sAIF, theta_Selectors_sAIF);
eta_Masks_sAIF = ismember(1:N_eta_sAIF, eta_Selectors_sAIF);
mu_Masks_sAIF = ismember(1:N_mu_sAIF, mu_Selectors_sAIF);

Means_Filtered_sAIF = Means_sAIF(delta_Masks_sAIF, theta_Masks_sAIF, eta_Masks_sAIF, mu_Masks_sAIF);
CVs_Filtered_sAIF = CVs_sAIF(delta_Masks_sAIF, theta_Masks_sAIF, eta_Masks_sAIF, mu_Masks_sAIF);

Means_Data_sAIF = reshape(Means_Filtered_sAIF, [], 1);
CVs_Data_sAIF = reshape(CVs_Filtered_sAIF, [], 1);
[delta_Grid_sAIF, theta_Grid_sAIF, eta_Grid_sAIF, mu_Grid_sAIF] = ndgrid(delta_vector_sAIF(delta_Masks_sAIF), theta_vector_sAIF(theta_Masks_sAIF), eta_vector_sAIF(eta_Masks_sAIF), mu_vector_sAIF(mu_Masks_sAIF));
delta_Data_sAIF = delta_Grid_sAIF(:);
theta_Data_sAIF = theta_Grid_sAIF(:);
eta_Data_sAIF = eta_Grid_sAIF(:);
mu_Data_sAIF = mu_Grid_sAIF(:);

%% Load and Process FP Data
Data_FP = load('FP_TwoSpecies_delta_theta.mat');
delta_vector_FP = Data_FP.delta_vector;
theta_vector_FP = Data_FP.theta_vector;
N_delta_FP = Data_FP.N_delta;
N_theta_FP = Data_FP.N_theta;
Means_FP = Data_FP.StationaryMean;
CVs_FP = sqrt(Data_FP.StationaryVariance) ./ Data_FP.StationaryMean;

% Filter Data
delta_Selectors_FP = 1:N_delta_FP;
theta_Selectors_FP = 1:N_theta_FP;

delta_Masks_FP = ismember(1:N_delta_FP, delta_Selectors_FP);
theta_Masks_FP = ismember(1:N_theta_FP, theta_Selectors_FP);

Means_Filtered_FP = Means_FP(delta_Masks_FP, theta_Masks_FP);
CVs_Filtered_FP = CVs_FP(delta_Masks_FP, theta_Masks_FP);

Means_Data_FP = reshape(Means_Filtered_FP, [], 1);
CVs_Data_FP = reshape(CVs_Filtered_FP, [], 1);

[delta_Grid_FP, theta_Grid_FP] = ndgrid(delta_vector_FP(delta_Masks_FP), theta_vector_FP(theta_Masks_FP));
delta_Data_FP = delta_Grid_FP(:);
theta_Data_FP = theta_Grid_FP(:);

%% Load and Process FP Data with high theta
Data_FP_Hightheta = load('FP_TwoSpecies_FixedHightheta.mat');
delta_vector_Hightheta = Data_FP_Hightheta.delta_vector;
N_delta_FP_Hightheta = Data_FP_Hightheta.N_delta;
Means_FP_Hightheta = Data_FP_Hightheta.StationaryMean;
CVs_FP_Hightheta = sqrt(Data_FP_Hightheta.StationaryVariance) ./ Data_FP_Hightheta.StationaryMean;

% Filter Data
delta_Selectors_Hightheta = 1:N_delta_FP_Hightheta;

delta_Masks_FP_Hightheta = ismember(1:N_delta_FP_Hightheta, delta_Selectors_Hightheta);

Means_Filtered_FP_Hightheta = Means_FP_Hightheta(delta_Masks_FP_Hightheta);
CVs_Filtered_FP_Hightheta = CVs_FP_Hightheta(delta_Masks_FP_Hightheta);

Means_Data_FP_Hightheta = reshape(Means_Filtered_FP_Hightheta, [], 1);
CVs_Data_FP_Hightheta = reshape(CVs_Filtered_FP_Hightheta, [], 1);

delta_Grid_FP_Hightheta = delta_vector_Hightheta(delta_Masks_FP_Hightheta);
delta_Data_Hightheta = delta_Grid_FP_Hightheta(:);

%% Figure Settings
SS = 4;
Figure_Width = 7.5 * SS;
Figure_Height = 3 * SS;
FontSize = 5 * SS;
FontSize_Small = 3 * SS;
FontSize_Large = 6 * SS;
LineWidth = 4;
LineWidth_Thin = 0.01;
MarkerSize = 25;
NominalColors = [55,126,184; ... % Blue
                 77,175,74; ... % Green
                 255,127,0; ... % Orange
                 228,26,28; ... % Red
                 163, 124, 182]/255; ... % Purple

%% Set Figure 1
Figure1_Name = 'sAIF_TwoSpecies_delta_theta_eta_mu';
Handle_Figure1 = figure;
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [0, 0, Figure_Width/2/1.18, Figure_Height];
    Handle_Figure1.PaperPositionMode = 'auto';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];

%% Set Axis 1
Handle_Axis1 = axes(Handle_Figure1);
    Handle_Axis1.Position = [0.15, 0.1, 0.79, 0.87];
    Handle_Axis1.Box = 'on';
    Handle_Axis1.FontSize = FontSize;
    hold(Handle_Axis1, 'on');
    grid(Handle_Axis1, 'on');
    Handle_Axis1.Layer = 'top';
    Handle_Axis1.XMinorGrid = 'off';
    Handle_Axis1.YMinorGrid = 'off';
    Handle_Axis1.XScale = 'Log';
    Handle_Axis1.YScale = 'Linear';
    Handle_Axis1.XLim = [0.001, 0.1];
    Handle_Axis1.YLim = [0, 20];
    Handle_Axis1.ZLim = [0.2, 1.2];
% sAIF
[delta_Indices_sAIF, theta_Indices_sAIF, eta_Indices_sAIF, mu_Indices_sAIF] = ndgrid(1:sum(delta_Masks_sAIF), 1:sum(theta_Masks_sAIF), 1:sum(eta_Masks_sAIF), 1:sum(mu_Masks_sAIF));

% Grid
N_x = 15; N_y = 30; N_z = 30;
x_Grid = logspace(-3, -1, N_x);
y_Grid = linspace(0, 20, N_y);
z_Grid = linspace(0.2, 1.2, N_z);
[x_Grid, y_Grid, z_Grid] = meshgrid(x_Grid, y_Grid, z_Grid);
Selected_x = zeros(numel(x_Grid),1);
Selected_y = zeros(numel(x_Grid),1);
Selected_z = zeros(numel(x_Grid),1);
for i = 1 : numel(x_Grid)
    Grid_Point = [x_Grid(i), y_Grid(i), z_Grid(i)];
    Distances = sqrt( (Grid_Point(1) - delta_Data_sAIF).^2 + (Grid_Point(2) - Means_Data_sAIF).^2 + (Grid_Point(3) - CVs_Data_sAIF).^2 );
    [~, Index] = min(Distances);
    Selected_x(i) = delta_Data_sAIF(Index);
    Selected_y(i) = Means_Data_sAIF(Index);
    Selected_z(i) = CVs_Data_sAIF(Index);
end
% scatter3(Handle_Axis1, Selected_x, Selected_y, Selected_z, MarkerSize, 'MarkerFaceColor', 'flat', 'MarkerEdgeColor', 'k', 'CData', NominalColors(1,:), 'Marker', 'o', 'LineWidth', LineWidth_Thin, 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 0.1);

scatter3(Handle_Axis1, delta_Data_sAIF, Means_Data_sAIF, CVs_Data_sAIF, MarkerSize, 'MarkerFaceColor', 'flat', 'MarkerEdgeColor', 'k', 'CData', NominalColors(1,:), 'Marker', 'o', 'LineWidth', LineWidth_Thin, 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 0.4);
Handle_Axis1.View = [80, 10];

% FP
if surface == true
    Handle_Surface_Up = surf(Handle_Axis1, delta_Grid_FP, Means_Filtered_FP, CVs_Filtered_FP+0.028); % Add small value so that the surface wouldn't appear perforated 
        Handle_Surface_Up.FaceColor = NominalColors(5,:);
        Handle_Surface_Up.FaceAlpha = 0.7;
        Handle_Surface_Up.EdgeColor = NominalColors(5,:);
        Handle_Surface_Up.EdgeAlpha = 1;
else
    scatter3(Handle_Axis1, delta_Data_FP, Means_Data_FP, CVs_Data_FP, MarkerSize*5, 'MarkerFaceColor', 'flat', 'MarkerEdgeColor', 'none', 'CData', NominalColors(5,:), 'Marker', 'o', 'LineWidth', LineWidth_Thin);
end

% FP Wide delta
% [~, x_Grid] = ndgrid(delta_vector_Hightheta(delta_Masks_FP_Hightheta)', delta_vector_FP(delta_Masks_FP)');
% Means_Filtered_FP_Hightheta = repmat(Means_Data_FP_Hightheta, 1, sum(delta_Masks_FP));
% CVs_Filtered_FP_Hightheta = repmat(CVs_Data_FP_Hightheta, 1, sum(delta_Masks_FP));
% Means_Data_FP_Hightheta = Means_Filtered_FP_Hightheta(:);
% CVs_Data_FP_Hightheta = CVs_Filtered_FP_Hightheta(:);
% 
% if surface == true
%     Handle_Surface_Down = surf(Handle_Axis1, x_Grid, Means_Filtered_FP_Hightheta, CVs_Filtered_FP_Hightheta);
%         Handle_Surface_Down.FaceColor = 0.5*[1,1,1];
%         Handle_Surface_Down.FaceAlpha = 0.7;
%         Handle_Surface_Down.EdgeColor = 0.5*[1, 1, 1];
%         Handle_Surface_Down.EdgeAlpha = 1;
% else
%     scatter3(Handle_Axis1, x_Grid(:), Means_Filtered_FP_Hightheta(:), CVs_Filtered_FP_Hightheta(:), MarkerSize*5, 'MarkerFaceColor', 'flat', 'MarkerEdgeColor', 'none', 'CData', cmap(1,:), 'Marker', 'o', 'LineWidth', LineWidth_Thin);
% end

%% Save Figures
if Save_Flag == 1  
    Handle_Figure1.Color = [1, 1, 1];
    % Handle_Figure1.Color = 'none';
    % set(Handle_Figure1, 'InvertHardCopy', 'off');
  	% print(Handle_Figure1, Figure1_Name, '-dpdf', '-vector','-bestfit');
    print(Handle_Figure1, Figure1_Name, '-dpng', '-r600');
    Handle_Figure1.Color = [1, 1, 1];
end