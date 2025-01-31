%% Clear Workspace
close all;
clear;
clc;

Save_Flag = 0;

%% Fixed Parameters
gamma_1 = 1;
K_S = 1;
mu = 10;
theta = 1;
u_bar = gamma_1 * mu / theta;

%% Cutoff Frequency and Nominal Integral Gains
omega_0_vector = [0.6, 1.5, 20];
K_I_vector = [0.108, 0.175, 0.243];

%% Transfer Functions
n = length(omega_0_vector);
LoopGain = cell(n, 1);
for i = 1:n
    omega_0 = omega_0_vector(i);
    a2 = omega_0+gamma_1;
    a1 = omega_0*gamma_1;
    LoopGain{i} = tf([0 0 0 K_S*omega_0],...
        [1 a2 a1 0]);
end

%% Root Locus
K_max_vector = omega_0_vector * u_bar / (4 * mu);
N_K = 5000;
K = [linspace(0, K_max_vector(1), N_K); ...
     linspace(0, K_max_vector(2), N_K); ...
     linspace(0, K_max_vector(3), N_K)];
R = zeros(3, length(K), n);
for i = 1 : n
    R(:, :, i) = rlocus(LoopGain{i}, K(i,:));
end

%% Breaking Points
s_b = (sqrt(omega_0_vector.^2 + gamma_1^2 - gamma_1*omega_0_vector) - (omega_0_vector + gamma_1)) / 3;

%% Nominal Poles
Poles_Nominal = zeros(3, length(K_I_vector));
for i = 1 : length(K_I_vector)
    omega_0 = omega_0_vector(i);
    K_I = K_I_vector(i);
    Poles_Nominal(:,i) = roots([1, gamma_1 + omega_0, omega_0*gamma_1, omega_0*K_I*K_S]);
end

%% Figure Settings
ScalingFactor = 1.25;
SS = 4;
Figure_Width = 8 * SS;
Figure_Height = 3 * SS;
FontSize = ScalingFactor*5 * SS;
FontSize_Small = ScalingFactor*3 * SS;
FontSize_Large = ScalingFactor*6 * SS;
LineWidth = ScalingFactor*0.65 * SS;
LineWidth_Thick = ScalingFactor*1 * SS;
LineWidth_Thin = ScalingFactor*0.5 * SS;
MarkerSize = ScalingFactor*5 * SS;
MarkerSize_Small = ScalingFactor*4 * SS;
IColor1 = [55,126,184]/255; % Blue
IColor2 = [77,175,74]/255; % Green
IColor3 = [228,26,28]/255; % Red

%% Set Figure 1
Figure1_Name = 'BD_FilteredI_RootLocus1';
Handle_Figure1 = figure();
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [2, 2, Figure_Width/2, Figure_Height*1.2];
    Handle_Figure1.PaperPositionMode = 'auto';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];
    
%% Set Figure 2
Figure2_Name = 'BD_FilteredI_RootLocus2';
Handle_Figure2 = figure();
    Handle_Figure2.PaperUnits = 'centimeters';
    Handle_Figure2.Units = 'centimeters';
    Handle_Figure2.Position = [2, 2, Figure_Width/2, Figure_Height*1.2];
    Handle_Figure2.PaperPositionMode = 'auto';
    Handle_Figure2.PaperSize = [Handle_Figure2.PaperPosition(3), Handle_Figure2.PaperPosition(4)];
    
%% Set Figure 3
Figure3_Name = 'BD_FilteredI_RootLocus3';
Handle_Figure3 = figure();
    Handle_Figure3.PaperUnits = 'centimeters';
    Handle_Figure3.Units = 'centimeters';
    Handle_Figure3.Position = [2, 2, Figure_Width/2, Figure_Height*1.2];
    Handle_Figure3.PaperPositionMode = 'auto';
    Handle_Figure3.PaperSize = [Handle_Figure3.PaperPosition(3), Handle_Figure3.PaperPosition(4)];

%% Axis 1
Handle_Axis1 = axes(Handle_Figure1);
    Handle_Axis1.Position = [0.04, 0.16, 0.78, 0.76];
    Handle_Axis1.Box = 'on';
    Handle_Axis1.BoxStyle = 'full';
    Handle_Axis1.LineWidth = LineWidth_Thin;
    Handle_Axis1.FontSize = FontSize_Small*1.25;
    hold(Handle_Axis1, 'on');
    grid(Handle_Axis1, 'on');
    Handle_Axis1.XMinorGrid = 'off';
    Handle_Axis1.YMinorGrid = 'off';
    Handle_Axis1.XLim = [-2*gamma_1, 0.025];
    Handle_Axis1.YLim = [-1, 1];
    Handle_Axis1.XLabel.String = 'Real Axis';
    Handle_Axis1.XLabel.FontSize = FontSize;
    Handle_Axis1.YLabel.String = 'Imaginary Axis';
    Handle_Axis1.YLabel.FontSize = FontSize;
    Handle_Axis1.XTick = [Handle_Axis1.XLim(1), -gamma_1, -gamma_1/2, 0];
    Handle_Axis1.XTickLabel = {num2str(Handle_Axis1.XLim(1)), '$-\gamma_1$', '$-\frac{\gamma_1}{2}$', '0'};
    Handle_Axis1.TickLabelInterpreter = 'latex';
    Handle_Axis1.Title.String = ['$\omega_0 = ', num2str(omega_0_vector(1)), '$'];
    Handle_Axis1.Title.Interpreter = 'latex';
    Handle_Axis1.Title.FontSize = FontSize;
    Handle_Axis1.YAxisLocation = 'Right';
    set(Handle_Axis1, 'ColorScale', 'linear');
    colormap(Handle_Axis1, 'turbo');
for i = 1 : size(R,1)
	scatter(Handle_Axis1, real(R(i,:,1)), imag(R(i,:,1)), 3*MarkerSize, K(1,:), 'filled');
end
plot(Handle_Axis1, 0 * ones(1,2), Handle_Axis1.YLim, 'Color', 'r', 'LineStyle', '-', 'LineWidth', LineWidth_Thin);
plot(Handle_Axis1, real(Poles_Nominal(:,1)), imag(Poles_Nominal(:,1)), 'Marker', 'o', 'MarkerSize', MarkerSize_Small, 'Color', 'k', 'LineStyle', 'none', 'MarkerFaceColor', IColor1);
plot(Handle_Axis1, s_b(1) * [1, 1], Handle_Axis1.YLim, 'LineStyle', '--', 'Color', 'k', 'LineWidth', LineWidth_Thin);

%% Axis 2
Handle_Axis2 = axes(Handle_Figure2);
    Handle_Axis2.Position = [0.04, 0.16, 0.78, 0.76];
    Handle_Axis2.Box = 'on';
    Handle_Axis2.BoxStyle = 'full';
    Handle_Axis2.LineWidth = LineWidth_Thin;
    Handle_Axis2.FontSize = FontSize_Small*1.25;
    hold(Handle_Axis2, 'on');
    grid(Handle_Axis2, 'on');
    Handle_Axis2.XMinorGrid = 'off';
    Handle_Axis2.YMinorGrid = 'off';
    Handle_Axis2.XLim = [-2*gamma_1, 0.025];
    Handle_Axis2.YLim = [-1, 1];
    Handle_Axis2.XLabel.String = 'Real Axis';
    Handle_Axis2.XLabel.FontSize = FontSize;
    Handle_Axis2.YLabel.String = 'Imaginary Axis';
    Handle_Axis2.YLabel.FontSize = FontSize;
    Handle_Axis2.XTick = [Handle_Axis1.XLim(1), -gamma_1, -gamma_1/2, 0];
    Handle_Axis2.XTickLabel = {num2str(Handle_Axis1.XLim(1)), '$-\gamma_1$', '$-\frac{\gamma_1}{2}$', '0'};
    Handle_Axis2.TickLabelInterpreter = 'latex';
    Handle_Axis2.Title.String = ['$\omega_0 = ', num2str(omega_0_vector(2)), '$'];
    Handle_Axis2.Title.Interpreter = 'latex';
    Handle_Axis2.Title.FontSize = FontSize;
    Handle_Axis2.YAxisLocation = 'Right';
    set(Handle_Axis2, 'ColorScale', 'linear');
    colormap(Handle_Axis2, 'turbo');
for i = 1 : size(R,1)
	scatter(Handle_Axis2, real(R(i,:,2)), imag(R(i,:,2)), 3*MarkerSize, K(2,:), 'filled');
end
plot(Handle_Axis2, 0 * ones(1,2), Handle_Axis2.YLim, 'Color', 'r', 'LineStyle', '-', 'LineWidth', LineWidth_Thin);
plot(Handle_Axis2, real(Poles_Nominal(:,2)), imag(Poles_Nominal(:,2)), 'Marker', 'o', 'MarkerSize', MarkerSize_Small, 'Color', 'k', 'LineStyle', 'none', 'MarkerFaceColor', IColor2);
plot(Handle_Axis2, s_b(2) * [1, 1], Handle_Axis1.YLim, 'LineStyle', '--', 'Color', 'k', 'LineWidth', LineWidth_Thin);

%% Axis 3
Handle_Axis3 = axes(Handle_Figure3);
    Handle_Axis3.Position = [0.04, 0.16, 0.78, 0.76];
    Handle_Axis3.Box = 'on';
    Handle_Axis3.BoxStyle = 'full';
    Handle_Axis3.LineWidth = LineWidth_Thin;
    Handle_Axis3.FontSize = FontSize_Small*1.25;
    hold(Handle_Axis3, 'on');
    grid(Handle_Axis3, 'on');
    Handle_Axis3.XMinorGrid = 'off';
    Handle_Axis3.YMinorGrid = 'off';
    Handle_Axis3.XLim = [-2*gamma_1, 0.025];
    Handle_Axis3.YLim = [-1, 1];
    Handle_Axis3.XLabel.String = 'Real Axis';
    Handle_Axis3.XLabel.FontSize = FontSize;
    Handle_Axis3.YLabel.String = 'Imaginary Axis';
    Handle_Axis3.YLabel.FontSize = FontSize;
    Handle_Axis3.XTick = [Handle_Axis1.XLim(1), -gamma_1, -gamma_1/2, 0];
    Handle_Axis3.XTickLabel = {num2str(Handle_Axis1.XLim(1)), '$-\gamma_1$', '$-\frac{\gamma_1}{2}$', '0'};
    Handle_Axis3.TickLabelInterpreter = 'latex';
    Handle_Axis3.Title.String = ['$\omega_0 = ', num2str(omega_0_vector(3)), '$'];
    Handle_Axis3.Title.Interpreter = 'latex';
    Handle_Axis3.Title.FontSize = FontSize;
    Handle_Axis3.YAxisLocation = 'Right';
    set(Handle_Axis3, 'ColorScale', 'linear');
    colormap(Handle_Axis3, 'turbo');
for i = 1 : size(R,1)
	scatter(Handle_Axis3, real(R(i,:,3)), imag(R(i,:,3)), 3*MarkerSize, K(3,:), 'filled');
end
plot(Handle_Axis3, 0 * ones(1,2), Handle_Axis3.YLim, 'Color', 'r', 'LineStyle', '-', 'LineWidth', LineWidth_Thin);
plot(Handle_Axis3, real(Poles_Nominal(:,3)), imag(Poles_Nominal(:,3)), 'Marker', 'o', 'MarkerSize', MarkerSize_Small, 'Color', 'k', 'LineStyle', 'none', 'MarkerFaceColor', IColor3);
plot(Handle_Axis3, s_b(3) * [1, 1], Handle_Axis1.YLim, 'LineStyle', '--', 'Color', 'k', 'LineWidth', LineWidth_Thin);

%% Save Figures
if Save_Flag == 1  
    Handle_Figure1.Color = 'none';
    set(Handle_Figure1, 'InvertHardCopy', 'off');
  	print(Handle_Figure1, Figure1_Name, '-dpdf', '-painters','-bestfit');
    Handle_Figure1.Color = [1, 1, 1];
    Handle_Figure2.Color = 'none';
    set(Handle_Figure2, 'InvertHardCopy', 'off');
    print(Handle_Figure2, Figure2_Name, '-dpdf', '-painters','-bestfit');
    Handle_Figure2.Color = [1, 1, 1];
    Handle_Figure3.Color = 'none';
    set(Handle_Figure3, 'InvertHardCopy', 'off');
    print(Handle_Figure3, Figure3_Name, '-dpdf', '-painters','-bestfit');
    Handle_Figure3.Color = [1, 1, 1];
end