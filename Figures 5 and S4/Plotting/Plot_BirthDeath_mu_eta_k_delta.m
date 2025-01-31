%% Clear Workspace
clear;
close all;
clc;
Save_Flag = 0;

%% Load Data
% Filtered-Proportional Controller
Data_FP = load('FP_BirthDeath_delta.mat');
Means_FP = Data_FP.StationaryMean;
CVs_FP = sqrt(Data_FP.StationaryVariance) ./ Data_FP.StationaryMean;
% rAIF Controller
Data_rAIF = load('rAIF_BirthDeath_mu_eta_k.mat');
Means_rAIF = Data_rAIF.StationaryMean;
CVs_rAIF = sqrt(Data_rAIF.StationaryVariance) ./ Data_rAIF.StationaryMean;
% sAIF Controller
Data_sAIF = load('sAIF_BirthDeath_mu_eta.mat');
Means_sAIF = Data_sAIF.StationaryMean;
CVs_sAIF = sqrt(Data_sAIF.StationaryVariance) ./ Data_sAIF.StationaryMean;
% Plant
gamma_1 = Data_FP.Parameters.gamma_1;

%% Analytical Formulas
Means_OL = 10.^linspace(-2, 2, 1e3);
CVs_OL = (1 ./ sqrt(Means_OL));

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
                 255,127,0; ... % Orange
                 228,26,28; ... % Red
                 163, 124, 182]/255; ... % Purple

%% Set Figure 1
Figure1_Name = 'BD_rAIF_Noise';
Handle_Figure1 = figure;
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [0, 0, Figure_Width/2/1.18, Figure_Height];
    Handle_Figure1.PaperPositionMode = 'auto';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];

%% Set Figure 2
Figure2_Name = 'BD_sAIF_Noise';
Handle_Figure2 = figure;
    Handle_Figure2.PaperUnits = 'centimeters';
    Handle_Figure2.Units = 'centimeters';
    Handle_Figure2.Position = [0, Figure_Height, Figure_Width/2, Figure_Height];
    Handle_Figure2.PaperPositionMode = 'auto';
    Handle_Figure2.PaperSize = [Handle_Figure2.PaperPosition(3), Handle_Figure2.PaperPosition(4)];

%% Set Axis 1: Topology I
Handle_Axis1 = axes(Handle_Figure1);
    Handle_Axis1.Position = [0.18, 0.16, 0.78, 0.76];
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
%     Handle_Axis1.XLabel.String = 'Stationary $E[X_2]$';
%     Handle_Axis1.YLabel.String = 'Stationary $CV[X_2]$';
    Handle_Axis1.YLabel.Interpreter = 'latex';
    Handle_Axis1.Title.String = 'Noise Amplification';
    Handle_Axis1.YLim = [0, 1];
    Handle_Axis1.XLim = [1, 10];
N_mu = size(Means_sAIF, 1);
N_eta = size(Means_sAIF, 2);
for i = 1:N_mu
    for j = 1:N_eta
        means = squeeze(Means_rAIF(i,j,:));
        CVs = squeeze(CVs_rAIF(i,j,:));
        Handle_rAIF = plot(Handle_Axis1, means, CVs, 'MarkerSize', MarkerSize, 'Marker', 'o', 'LineStyle', 'none', 'MarkerFaceColor', NominalColors(4,:), 'Color', 'k', 'LineWidth', LineWidth_Thin);
    end
end
Handle_OL = plot(Handle_Axis1, Means_OL, CVs_OL, 'k-', 'LineWidth', LineWidth);
legend(Handle_Axis1, [Handle_rAIF, Handle_OL], {'rAIF', 'Open Loop'});

%% Set Axis 2: Topology II
Handle_Axis2 = axes(Handle_Figure2);
    Handle_Axis2.Position = [0.15, 0.16, 0.8, 0.76];
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
%     Handle_Axis2.XLabel.String = 'Stationary $E[X_2]$';
%     Handle_Axis2.YLabel.String = 'Stationary $CV[X_2]$';
    Handle_Axis2.YLabel.Interpreter = 'latex';
    Handle_Axis2.Title.String = 'Noise Reduction';
    Handle_Axis2.YLim = [0, 1];
    Handle_Axis2.XLim = [1, 10];
N_eta = size(Means_sAIF, 2);
cmap = colormap(turbo(N_eta));
for i = 1:N_eta
    means = Means_sAIF(:,i);
    CVs = CVs_sAIF(:,i);
    if i == 3
        Handle_PI1 = plot(Handle_Axis2, means, CVs, 'MarkerSize', MarkerSize, 'Marker', 'o', 'LineStyle', 'none', 'MarkerFaceColor', cmap(i,:), 'Color', 'k', 'LineWidth', LineWidth_Thin);
    elseif i == 5
        Handle_PI2 = plot(Handle_Axis2, means, CVs, 'MarkerSize', MarkerSize, 'Marker', 'o', 'LineStyle', 'none', 'MarkerFaceColor', cmap(i,:), 'Color', 'k', 'LineWidth', LineWidth_Thin);
    else
        plot(Handle_Axis2, means, CVs, 'MarkerSize', MarkerSize, 'Marker', 'o', 'LineStyle', 'none', 'MarkerFaceColor', cmap(i,:), 'Color', 'k', 'LineWidth', LineWidth_Thin);
    end
end
clim([Data_sAIF.eta_vector(1), Data_sAIF.eta_vector(end)]);
Handle_Colorbar = colorbar();
    % Handle_Colorbar.Label.String = '$\eta$';
    % Handle_Colorbar.Label.Interpreter = 'latex';
    % Handle_Colorbar.Label.Position = [1.5, 50, 0];
    % Handle_Colorbar.Label.FontSize = FontSize*1.5;
    Handle_Colorbar.Ruler.Scale = 'log';
    Handle_Colorbar.Ticks = [Data_sAIF.eta_vector(1), Data_sAIF.eta_vector(end)];
Handle_OL = plot(Handle_Axis2, Means_OL, CVs_OL, 'k-', 'LineWidth', LineWidth);
Handle_P = plot(Handle_Axis2, reshape(Means_FP,1,[]), reshape(CVs_FP,1,[]), 'MarkerSize', MarkerSize*1.5, 'Marker', 'o', 'LineStyle', 'none', 'MarkerFaceColor', NominalColors(5,:), 'Color', 'k', 'LineWidth', LineWidth_Thin);
[h, icons] = legend(Handle_Axis2, [Handle_PI1, Handle_OL, Handle_P], {'sAIF', 'Open Loop', 'Filtered P'});

%% Save Figures
if Save_Flag == 1  
    Handle_Figure1.Color = 'none';
    set(Handle_Figure1, 'InvertHardCopy', 'off');
  	print(Handle_Figure1, Figure1_Name, '-dpdf', '-vector','-bestfit');
    Handle_Figure1.Color = [1, 1, 1];
    Handle_Figure2.Color = 'none';
    set(Handle_Figure2, 'InvertHardCopy', 'off');
    print(Handle_Figure2, Figure2_Name, '-dpdf', '-vector','-bestfit');
    Handle_Figure2.Color = [1, 1, 1];
    Handle_Figure3.Color = 'none';
end