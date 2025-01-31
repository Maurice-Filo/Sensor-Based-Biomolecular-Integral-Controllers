%% Clear Workspace
clear;
close all;
clc;
Save_Flag = 0;
ColorParameter = 'eta';

%% Load and Process sAIF Data
Data_sAIF_mu_1 = load('sAIF_BirthDeath_delta_theta_eta_mu1.mat');
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
    Data_sAIF_mu_i = load(['sAIF_BirthDeath_delta_theta_eta_mu', num2str(i), '.mat']);
    Means_sAIF(:,:,:,i) = Data_sAIF_mu_i.StationaryMean;
    CVs_sAIF(:,:,:,i) = sqrt(Data_sAIF_mu_i.StationaryVariance) ./ Data_sAIF_mu_i.StationaryMean;
end

% Filter Data
delta_Selectors_sAIF = 9;%1:N_delta_sAIF;
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
[delta_Data_sAIF, theta_Data_sAIF, eta_Data_sAIF, mu_Data_sAIF] = ndgrid(delta_vector_sAIF(delta_Masks_sAIF), theta_vector_sAIF(theta_Masks_sAIF), eta_vector_sAIF(eta_Masks_sAIF), mu_vector_sAIF(mu_Masks_sAIF));
delta_Data_sAIF = delta_Data_sAIF(:);
theta_Data_sAIF = theta_Data_sAIF(:);
eta_Data_sAIF = eta_Data_sAIF(:);
mu_Data_sAIF = mu_Data_sAIF(:);

%% Load and Process FP Data
Data_FP = load('FP_BirthDeath_delta_theta.mat');
delta_vector_FP = Data_FP.delta_vector;
theta_vector_FP = Data_FP.theta_vector;
N_delta_FP = Data_FP.N_delta;
N_theta_FP = Data_FP.N_theta;
Means_FP = Data_FP.StationaryMean;
CVs_FP = sqrt(Data_FP.StationaryVariance) ./ Data_FP.StationaryMean;

% Filter Data
delta_Selectors_FP = 9; %1:N_delta_FP;
theta_Selectors_FP = 1:N_theta_FP;

delta_Masks_FP = ismember(1:N_delta_FP, delta_Selectors_FP);
theta_Masks_FP = ismember(1:N_theta_FP, theta_Selectors_FP);

Means_Filtered_FP = Means_FP(delta_Masks_FP, theta_Masks_FP);
CVs_Filtered_FP = CVs_FP(delta_Masks_FP, theta_Masks_FP);

Means_Data_FP = reshape(Means_Filtered_FP, [], 1);
CVs_Data_FP = reshape(CVs_Filtered_FP, [], 1);

[delta_Data_FP, theta_Data_FP] = ndgrid(delta_vector_FP(delta_Masks_FP), theta_vector_FP(theta_Masks_FP));
delta_Data_FP = delta_Data_FP(:);
theta_Data_FP = theta_Data_FP(:);

%% Load and Process FP Data with high theta
Data_FP_Hightheta = load('FP_BirthDeath_FixedHightheta.mat');
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
Figure_Width = 6 * SS;
Figure_Height = 1.35 * SS;
FontSize = 5 * SS;
FontSize_Small = 3 * SS;
FontSize_Large = 6 * SS;
LineWidth = 4;
LineWidth_Thin = 0.01;
MarkerSize = 20;
NominalColors = [55,126,184; ... % Blue
                 77,175,74; ... % Green
                 255,127,0; ... % Orange
                 228,26,28; ... % Red
                 163, 124, 182]/255; ... % Purple

%% Set Figure 1
Figure1_Name = 'sAIF_BirthDeath_deltaSlice_theta_eta_mu';
Handle_Figure1 = figure;
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [0, 0, Figure_Width/2/1.18, Figure_Height];
    Handle_Figure1.PaperPositionMode = 'auto';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];

%% Set Axis 1
Handle_Axis1 = axes(Handle_Figure1);
    Handle_Axis1.Position = [0.15, 0.17, 0.8, 0.77];
    Handle_Axis1.Box = 'on';
    Handle_Axis1.FontSize = FontSize;
    hold(Handle_Axis1, 'on');
    grid(Handle_Axis1, 'on');
    Handle_Axis1.Layer = 'top';
    Handle_Axis1.XMinorGrid = 'off';
    Handle_Axis1.YMinorGrid = 'off';
    Handle_Axis1.XScale = 'linear';
    Handle_Axis1.YScale = 'linear';
    Handle_Axis1.XLim = [0, 20];
    Handle_Axis1.YLim = [0.2, 1.2];
    Handle_Axis1.YTick = 0.2:0.2:1.2;
% sAIF
[delta_Indices_sAIF, theta_Indices_sAIF, eta_Indices_sAIF, mu_Indices_sAIF] = ndgrid(1:sum(delta_Masks_sAIF), 1:sum(theta_Masks_sAIF), 1:sum(eta_Masks_sAIF), 1:sum(mu_Masks_sAIF));
% switch ColorParameter
%     case 'delta'
%         cmap = colormap(turbo(sum(delta_Masks_sAIF)));
%         ColorIndices = delta_Indices_sAIF(:);
%         Colors = cmap(ColorIndices, :);
%         Range = [min(delta_vector_sAIF(delta_Selectors_sAIF)), max(delta_vector_sAIF(delta_Selectors_sAIF))];
%     case 'theta'
%         cmap = colormap(turbo(sum(theta_Masks_sAIF)));
%         ColorIndices = theta_Indices_sAIF(:);
%         Colors = cmap(ColorIndices, :);
%         Range = [min(theta_vector_sAIF(theta_Selectors_sAIF)), max(theta_vector_sAIF(theta_Selectors_sAIF))];
%     case 'eta'
%         cmap = colormap(turbo(sum(eta_Masks_sAIF)));
%         ColorIndices = eta_Indices_sAIF(:);
%         Colors = cmap(ColorIndices, :);
%         Range = [min(eta_vector_sAIF(eta_Selectors_sAIF)), max(eta_vector_sAIF(eta_Selectors_sAIF))];
%     case 'mu'
%         cmap = colormap(turbo(sum(mu_Masks_sAIF)));
%         ColorIndices = mu_Indices_sAIF(:);
%         Colors = cmap(ColorIndices, :);
%         Range = [min(mu_vector_sAIF(mu_Selectors_sAIF)), max(mu_vector_sAIF(mu_Selectors_sAIF))];
% end
scatter(Handle_Axis1, Means_Data_sAIF, CVs_Data_sAIF, MarkerSize, 'MarkerFaceColor', 'flat', 'CData', NominalColors(1,:), 'Marker', 'o', 'MarkerEdgeColor', 'k', 'LineWidth', LineWidth_Thin, 'MarkerEdgeAlpha', 0.4);
% clim(Handle_Axis1, Range);
% Handle_Colorbar = colorbar();
%     Handle_Colorbar.Ruler.Scale = 'log';
%     Handle_Colorbar.Ticks = Range;

% FP
plot(Handle_Axis1, Means_Data_FP, CVs_Data_FP, 'Color', NominalColors(5,:), 'LineWidth', LineWidth*1.5);

% FP Wide delta
% plot(Handle_Axis1, Means_Data_FP_Hightheta, CVs_Data_FP_Hightheta, 'Color', 0.5*[1, 1, 1], 'LineWidth', LineWidth*1.5);

%% Save Figures
if Save_Flag == 1  
    % Handle_Figure1.Color = 'none';
    % set(Handle_Figure1, 'InvertHardCopy', 'off');
  	% print(Handle_Figure1, Figure1_Name, '-dpdf', '-vector','-bestfit');
    print(Handle_Figure1, Figure1_Name, '-dpng', '-r600');
    Handle_Figure1.Color = [1, 1, 1];
end

%% Enable data cursor mode
dcm = datacursormode(gcf);
parameter_strings1 = {'\delta', '\theta', '\eta', '\mu'}; 
x1 = Means_Data_sAIF;
y1 = CVs_Data_sAIF;
W1 = [delta_Data_sAIF, theta_Data_sAIF, eta_Data_sAIF, mu_Data_sAIF];

parameter_strings2 = {'\delta', '\theta'}; 
x2 = Means_Data_FP;
y2 = CVs_Data_FP;
W2 = [delta_Data_FP, theta_Data_FP];

set(dcm, 'UpdateFcn', @(src, event) DisplayValues2(src, event, x1, y1, W1, parameter_strings1, x2, y2, W2, parameter_strings2));
