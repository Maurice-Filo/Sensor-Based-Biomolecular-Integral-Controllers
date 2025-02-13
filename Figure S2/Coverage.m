%% Clear Workspace
close all;
clear;
clc

%% Parameters
u_bar = 8; 
mu = 1;

%% Cooperativity for repression
n_1 = 1;
n_2 = 2;

%% Ranges
cubic_length = 10*n_1;
K_I = linspace(0, cubic_length, 10);
K_P = linspace(0, 2*cubic_length, 10);
omega_c = linspace(0, cubic_length, 10);

%% Plot Specs
Colors = lines(10);

% Color_whole = 'b';
% Alpha_whole = 0.1;
% EdgeColor_whole = 'b';

Color_deg = [207, 92, 54] / 255;
Alpha_deg = 0.1;
EdgeColor_deg = [207, 92, 54] / 255;

Color_rep1 = [81, 137, 188] / 255;
Alpha_rep1 = 0.2;
EdgeColor_rep1 = [81, 137, 188] / 255;

Color_rep = [78, 146, 103] / 255;
Alpha_rep = 0.3;
EdgeColor_rep = [78, 146, 103] / 255;

Cut = 0.01;
LineWidth = 3;
LineWidth_Thin = 1;

%% Figure Settings
ScalingFactor = 1;
SS = 4;
Figure_Width = 10 * SS;
Figure_Height = 5 * SS;
FontSize = ScalingFactor*8 * SS;
FontSize_Small = ScalingFactor*3 * SS;
FontSize_Large = ScalingFactor*6 * SS;

%% Set Figure 1
Figure1_Name = 'Coverages';
Handle_Figure1 = figure();
    Handle_Figure1.Color = [1, 1, 1];
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [SS, SS, Figure_Width/2, Figure_Height];
    Handle_Figure1.PaperPositionMode = 'auto';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];

%% Axis for Coverages
Handle_Axis1 = axes(Handle_Figure1);
    Handle_Axis1.Box = 'on';
    Handle_Axis1.BoxStyle = 'full';
    Handle_Axis1.LineWidth = LineWidth_Thin;
    Handle_Axis1.FontSize = FontSize;
    hold(Handle_Axis1, 'on');
    grid(Handle_Axis1, 'on');
    Handle_Axis1.XMinorGrid = 'off';
    Handle_Axis1.YMinorGrid = 'off';
    Handle_Axis1.XLabel.String = '$K_I$';
    Handle_Axis1.XLabel.Position = [5.826495203284537,22.347880863540837,-1.222694313458248];
    Handle_Axis1.YLabel.String = '$K_P$';
    Handle_Axis1.ZLabel.String = '$\omega_0$';
    Handle_Axis1.XLabel.Interpreter = 'latex';
    Handle_Axis1.YLabel.Interpreter = 'latex';
    Handle_Axis1.ZLabel.Interpreter = 'latex';
    Handle_Axis1.ZScale = 'linear';
    Handle_Axis1.ColorScale = 'linear';
    Handle_Axis1.View = [130,13];
    Handle_Axis1.TickLabelInterpreter = 'latex';
    Handle_Axis1.XTick = 0;
    Handle_Axis1.YTick = [0, n_1*u_bar/mu, n_2*u_bar/mu];
    Handle_Axis1.ZTick = 0;
    Handle_Axis1.XTickLabel = {'0'};
    Handle_Axis1.YTickLabel = {'0', '$\frac{\bar u}{\mu}$', '$\frac{2\bar u}{\mu}$'};
    Handle_Axis1.XLim = [0, 10+3*Cut];
    Handle_Axis1.YLim = [0, 20];
    Handle_Axis1.ZLim = [0, 10];
    Handle_Axis1.View = [109.024,   43.5919];

% Whole space
% Handle_Fill_whole_1 = fill3(Handle_Axis1, [K_I(1),         K_I(1),         K_I(end),       K_I(end),       K_I(1)], ...
%                                         [K_P(1),         K_P(end),       K_P(end),       K_P(1),         K_P(1)], ...
%                                         [omega_c(1),     omega_c(1),     omega_c(1),     omega_c(1),     omega_c(1)], Color_whole, 'FaceAlpha', Alpha_whole, 'EdgeColor', EdgeColor_whole, 'LineWidth', LineWidth);
% Handle_Fill_whole_2 = fill3(Handle_Axis1, [K_I(1),         K_I(1),         K_I(end),       K_I(end),       K_I(1)], ...
%                                         [K_P(1),         K_P(1),         K_P(1),         K_P(1),         K_P(1)], ...
%                                         [omega_c(1),     omega_c(end),   omega_c(end),   omega_c(1),     omega_c(1)], Color_whole, 'FaceAlpha', Alpha_whole, 'EdgeColor', EdgeColor_whole, 'LineWidth', LineWidth);
% Handle_Fill_whole_3 = fill3(Handle_Axis1, [K_I(1),         K_I(1),         K_I(end),       K_I(end),       K_I(1)], ...
%                                         [K_P(1),         K_P(end),       K_P(end),       K_P(1),         K_P(1)], ...
%                                         [omega_c(end),   omega_c(end),   omega_c(end),   omega_c(end),   omega_c(end)], Color_whole, 'FaceAlpha', Alpha_whole, 'EdgeColor', EdgeColor_whole, 'LineWidth', LineWidth);   
% Handle_Fill_whole_4 = fill3(Handle_Axis1, [K_I(1),         K_I(end),       K_I(end),       K_I(1),         K_I(1)], ...
%                                         [K_P(end),       K_P(end),       K_P(end),       K_P(end),       K_P(end)], ...
%                                         [omega_c(1),     omega_c(1),     omega_c(end),   omega_c(end),   omega_c(1)], Color_whole, 'FaceAlpha', Alpha_whole, 'EdgeColor', EdgeColor_whole, 'LineWidth', LineWidth);
% Handle_Fill_whole_5 = fill3(Handle_Axis1, [K_I(end),       K_I(end),       K_I(end),       K_I(end),       K_I(end)], ...
%                                         [K_P(1),         K_P(end),       K_P(end),       K_P(1),         K_P(1)], ...
%                                         [omega_c(1),     omega_c(1),     omega_c(end),   omega_c(end),   omega_c(1)], Color_whole, 'FaceAlpha', Alpha_whole, 'EdgeColor', EdgeColor_whole, 'LineWidth', LineWidth);
% Handle_Fill_whole_6 = fill3(Handle_Axis1, [K_I(1),         K_I(1),         K_I(1),         K_I(1),         K_I(1)], ...
%                                         [K_P(1),         K_P(end),       K_P(end),       K_P(1),         K_P(1)], ...
%                                         [omega_c(1),     omega_c(1),     omega_c(end),   omega_c(end),   omega_c(1)], Color_whole, 'FaceAlpha', Alpha_whole, 'EdgeColor', EdgeColor_whole, 'LineWidth', LineWidth);
  
% antithetic Filtered-PI Coverage Degradationf
K_P_trim_deg = linspace(K_I(end)/omega_c(end), cubic_length, 50);
omega_c_wrt_K_P = K_I(end) ./ K_P_trim_deg;

[y, z] = meshgrid(K_P, omega_c);
x = y.*z;
Handle_Surf_1 = surf(Handle_Axis1, x, y, z, 'FaceAlpha', Alpha_deg, 'FaceColor', Color_deg, 'EdgeColor', [0.5, 0.5, 0.5], 'LineWidth', 0.01);
% colormap hot

Handle_Fill_deg_1 = fill3(Handle_Axis1, [K_I(1),       K_I(end),           K_I(end),       K_I(1),         K_I(1)],...
                                         [K_P(end),     K_P(end),           K_P(end),       K_P(end),       K_P(end)],...
                                         [omega_c(1),   K_I(end)/K_P(end),  omega_c(end),   omega_c(end),   omega_c(1)],...
                                         Color_deg, 'FaceAlpha', Alpha_deg, 'EdgeColor', EdgeColor_deg, 'LineWidth', 2*LineWidth);
Handle_Fill_deg_2 = fill3(Handle_Axis1, [K_I(1),       K_I(1),         K_I(end),               K_I(end)],...
                                         [K_P(end),     K_P(1),         K_I(end)/omega_c(end),  K_P(end)],...
                                         [omega_c(end), omega_c(end),   omega_c(end),           omega_c(end)],...
                                         Color_deg, 'FaceAlpha', Alpha_deg, 'EdgeColor', EdgeColor_deg, 'LineWidth', 2*LineWidth);
Handle_Fill_deg_3 = fill3(Handle_Axis1, [K_I(end),             K_I(end),       K_I(end),               repmat(K_I(end), 1, length(K_P_trim_deg)),  K_I(end)],...
                                         [K_P(end),             K_P(end),       K_I(end)/omega_c(end),  K_P_trim_deg,                               K_P(end)],...
                                         [K_I(end)/K_P(end),    omega_c(end),   omega_c(end),           omega_c_wrt_K_P,                        K_I(end)/K_P(end)],...
                                         Color_deg, 'FaceAlpha', Alpha_deg, 'EdgeColor', EdgeColor_deg, 'LineWidth', 2*LineWidth);
Handle_Fill_deg_4 = fill3(Handle_Axis1, [K_I(1),       K_I(1),     K_I(1),         K_I(1),         K_I(1)],...
                                         [K_P(1),       K_P(end),   K_P(end),       K_P(1),         K_P(1)],...
                                         [omega_c(1),   omega_c(1), omega_c(end),   omega_c(end),   omega_c(1)],...
                                         Color_deg, 'FaceAlpha', Alpha_deg, 'EdgeColor', EdgeColor_deg, 'LineWidth', 2*LineWidth);
clear x y z

% antithetic Filtered-PI Coverage Repression n_2
a = n_2*u_bar/mu;
qua_sol = [(a-sqrt(a^2-4*a*K_I(end)/omega_c(end)))/2, (a+sqrt(a^2-4*a*K_I(end)/omega_c(end)))/2];
K_P_trim_quar = linspace(qua_sol(1), qua_sol(2), 50);
omega_c_wrt_K_P_trim = K_I(end)./(K_P_trim_quar.*(1-K_P_trim_quar./a));

% For better visualization
K_I(end) = K_I(end) + Cut;

K_P_trim_rep = linspace(0, a, 50);
[y, z] = meshgrid(K_P_trim_rep, omega_c);
x = z.*y.*(1-y/a);
Handle_Surf_21 = surf(Handle_Axis1, x, y, z, 'FaceAlpha', Alpha_rep, 'FaceColor', Color_rep1, 'EdgeColor', [0.5, 0.5, 0.5], 'LineWidth', 0.01);
% colormap summer

Handle_Fill_rep_11 = fill3(Handle_Axis1,  [K_I(1),       K_I(end),      K_I(end),      K_I(1),         K_I(1)],...
                                         [K_P(1),       qua_sol(1),    qua_sol(2),    a,              K_P(1)],...
                                         [omega_c(end), omega_c(end),  omega_c(end),  omega_c(end),   omega_c(end)],...
                                         Color_rep1, 'FaceAlpha', Alpha_rep1, 'EdgeColor', EdgeColor_rep1, 'LineWidth', 2*LineWidth);
Handle_Fill_rep_21 = fill3(Handle_Axis1,  [K_I(1),       K_I(1),         K_I(1),     K_I(1),     K_I(1)],...
                                         [K_P(1),       a,              a,          K_P(1),     K_P(1)],...
                                         [omega_c(end), omega_c(end),   omega_c(1), omega_c(1), omega_c(end)],...
                                         Color_rep1, 'FaceAlpha', Alpha_rep1, 'EdgeColor', EdgeColor_rep1, 'LineWidth', 2*LineWidth);
Handle_Fill_rep_31 = fill3(Handle_Axis1,  [K_I(end),         repmat(K_I(end), 1, length(K_P_trim_quar)),  K_I(end),       K_I(end)],...
                                         [qua_sol(1),       K_P_trim_quar,                               qua_sol(2),     qua_sol(1)],...
                                         [omega_c(end),     omega_c_wrt_K_P_trim,                       omega_c(end),   omega_c(end)],...
                                         Color_rep1, 'FaceAlpha', Alpha_rep1, 'EdgeColor', EdgeColor_rep1, 'LineWidth', 2*LineWidth);
Handle_Fill_rep_41 = fill3(Handle_Axis1, [K_I(1),       K_I(1),     K_I(1),         K_I(1),         K_I(1)],...
                                         [K_P(1),            a,          a,       K_P(1),         K_P(1)],...
                                         [omega_c(1),   omega_c(1), omega_c(end),   omega_c(end),   omega_c(1)],...
                                         Color_rep1, 'FaceAlpha', Alpha_rep1, 'EdgeColor', EdgeColor_rep1, 'LineWidth', 2*LineWidth);
clear x y z

% antithetic Filtered-PI Coverage Repression n_1
K_I(end) = K_I(end) + 2 * Cut;

a = n_1*u_bar/mu;
qua_sol = [(a-sqrt(a^2-4*a*K_I(end)/omega_c(end)))/2, (a+sqrt(a^2-4*a*K_I(end)/omega_c(end)))/2];
K_P_trim_quar = linspace(qua_sol(1), qua_sol(2), 50);
omega_c_wrt_K_P_trim = K_I(end)./(K_P_trim_quar.*(1-K_P_trim_quar./a));

% For better visualization
K_P_trim_rep = linspace(0, a, 50);
[y, z] = meshgrid(K_P_trim_rep, omega_c);
x = z.*y.*(1-y/a);
Handle_Surf_2 = surf(Handle_Axis1, x, y, z, 'FaceAlpha', Alpha_rep, 'FaceColor', Color_rep, 'EdgeColor', [0.5, 0.5, 0.5], 'LineWidth', 0.01);
% colormap summer

Handle_Fill_rep_1 = fill3(Handle_Axis1,  [K_I(1),       K_I(end),      K_I(end),      K_I(1),         K_I(1)],...
                                         [K_P(1),       qua_sol(1),    qua_sol(2),    a,              K_P(1)],...
                                         [omega_c(end), omega_c(end),  omega_c(end),  omega_c(end),   omega_c(end)],...
                                         Color_rep, 'FaceAlpha', Alpha_rep, 'EdgeColor', EdgeColor_rep, 'LineWidth', 2*LineWidth);
Handle_Fill_rep_2 = fill3(Handle_Axis1,  [K_I(1),       K_I(1),         K_I(1),     K_I(1),     K_I(1)],...
                                         [K_P(1),       a,              a,          K_P(1),     K_P(1)],...
                                         [omega_c(end), omega_c(end),   omega_c(1), omega_c(1), omega_c(end)],...
                                         Color_rep, 'FaceAlpha', Alpha_rep, 'EdgeColor', EdgeColor_rep, 'LineWidth', 2*LineWidth);
Handle_Fill_rep_3 = fill3(Handle_Axis1,  [K_I(end),         repmat(K_I(end), 1, length(K_P_trim_quar)),  K_I(end),       K_I(end)],...
                                         [qua_sol(1),       K_P_trim_quar,                               qua_sol(2),     qua_sol(1)],...
                                         [omega_c(end),     omega_c_wrt_K_P_trim,                       omega_c(end),   omega_c(end)],...
                                         Color_rep, 'FaceAlpha', Alpha_rep, 'EdgeColor', EdgeColor_rep, 'LineWidth', 2*LineWidth);
Handle_Fill_rep_4 = fill3(Handle_Axis1, [K_I(1),       K_I(1),     K_I(1),         K_I(1),         K_I(1)],...
                                         [K_P(1),           a,          a,       K_P(1),         K_P(1)],...
                                         [omega_c(1),   omega_c(1), omega_c(end),   omega_c(end),   omega_c(1)],...
                                         Color_rep, 'FaceAlpha', Alpha_rep, 'EdgeColor', EdgeColor_rep, 'LineWidth', 2*LineWidth);

                                     
%% Annotations
figure1 = Handle_Figure1;
annotation(figure1,'textarrow',[0.159833333333333 0.258611111111112],...
    [0.904194444444444 0.40325],...
    'String',{'$\omega_0 = \frac{K_I}{K_P(1-\frac{\mu K_P}{\bar{u}})}$'},...
    'Interpreter','latex',...
    'FontSize',28,...
    'FontName','Helvetica Neue');
annotation(figure1,'textarrow',[0.627263888888889 0.482625],...
    [0.0892777777777777 0.297416666666667],...
    'String',{'$\omega_0 = \frac{K_I}{K_P(1-\frac{\mu K_P}{2\bar{u}})}$'},...
    'Interpreter','latex',...
    'FontSize',28,...
    'FontName','Helvetica Neue');
annotation(figure1,'textarrow',[0.198638888888889 0.297416666666667],...
    [0.113972222222222 0.278013888888889],...
    'String',{'$\omega_0 = \frac{K_I}{K_P}$'},...
    'Interpreter','latex',...
    'FontSize',28,...
    'FontName','Helvetica Neue');

annotation(figure1,'textarrow',[0.4985 0.366208333333333],...
    [0.935944444444444 0.690763888888889],...
    'String',{'Repression $\mathcal{S}^1_r$'},...
    'Interpreter','latex',...
    'FontSize',20,...
    'FontName','Helvetica Neue');
annotation(figure1,'textarrow',[0.756027777777778 0.641375],...
    [0.888319444444444 0.713694444444444],...
    'String',{'Repression $\mathcal{S}^2_r$'},...
    'Interpreter','latex',...
    'FontSize',20,...
    'FontName','Helvetica Neue');

annotation(figure1,'textarrow',[0.793069444444444 0.743680555555556],...
    [0.278013888888889 0.375027777777778],...
    'String',{'Degradation $\mathcal{S}_d$'},...
    'Interpreter','latex',...
    'FontSize',20,...
    'FontName','Helvetica Neue');


%% Save Figures
Save_Flag = 1;
if Save_Flag == 1  
    set(Handle_Figure1, 'DefaultFigureRendererMode', 'manual');
    Handle_Figure1.Color = 'none';
    set(Handle_Figure1, 'InvertHardCopy', 'off');
    print(Handle_Figure1, Figure1_Name, '-dpdf', '-vector');
    Handle_Figure1.Color = [1, 1, 1];
end