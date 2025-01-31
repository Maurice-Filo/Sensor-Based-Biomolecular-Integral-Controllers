%% Clear Workspace
close all;
clear;
clc;

Save_Flag = 0;

%% Cuttoff Frequencies & Integral Gain
omega_0_vector = [0.6, 1.5, 20];
K_I_vector = [0.108, 0.175, 0.243];

%% Plant Parameters
Parameters.gamma_1 = 1;

%% Fixed Controller Parameters
Parameters.mu = 5;
Parameters.theta = 1;

%% Setpoint and Disturbance
r = 5;
Factor_Disturbance = 2;
r_Disturbed = r * Factor_Disturbance;

%% Functions
StoichiometryMatrix = StoichiometryMatrix_FilteredI_BD();
PropensityFunction = @PropensityFunction_FilteredI_BD;
SupportingInput = @SupportingInput_BD;
FixedPoint = @I_FixedPoint;
Params2Gains = @I_Params2Gains;
Gains2Params = @I_Gains2Params;
OutputIndex = 1;

%% Simulation Settings
N_Simulation = 1000;
Solver = 'ODE23s';
t_1 = 100;
t_f = 120;

%% Fixed Gains
Gains.K_S = Parameters.theta;

%% Simulations
% AIF Controller 1
Gains.omega_0 = omega_0_vector(1);
Gains.K_I = K_I_vector(1);
Parameters_Disturbed = Parameters;
Parameters_Disturbed.mu = Parameters.mu * Factor_Disturbance;
[Parameters_Disturbed, ~] = Gains2Params(Gains, Parameters_Disturbed, r_Disturbed, SupportingInput);
Parameters = Parameters_Disturbed; Parameters.mu = Parameters_Disturbed.mu / Factor_Disturbance;

IC = FixedPoint(Parameters, SupportingInput);
[time_vector_1, X_1] = DSA(StoichiometryMatrix, PropensityFunction, Parameters, IC, t_1, N_Simulation, Solver);
[time_vector_2, X_2] = DSA(StoichiometryMatrix, PropensityFunction, Parameters_Disturbed, X_1(:,end), t_f-t_1, N_Simulation, Solver);
X = [X_1, X_2(:,2:end)];
time_vector = [time_vector_1, t_1 + time_vector_2(2:end)];
x_AIF_1 = X(OutputIndex, :);

% AIF Controller 2
Gains.omega_0 = omega_0_vector(2);
Gains.K_I = K_I_vector(2);
Parameters_Disturbed = Parameters;
Parameters_Disturbed.mu = Parameters.mu * Factor_Disturbance;
[Parameters_Disturbed, ~] = Gains2Params(Gains, Parameters_Disturbed, r_Disturbed, SupportingInput);
Parameters = Parameters_Disturbed; Parameters.mu = Parameters_Disturbed.mu / Factor_Disturbance;

IC = FixedPoint(Parameters, SupportingInput);
[time_vector_1, X_1] = DSA(StoichiometryMatrix, PropensityFunction, Parameters, IC, t_1, N_Simulation, Solver);
[time_vector_2, X_2] = DSA(StoichiometryMatrix, PropensityFunction, Parameters_Disturbed, X_1(:,end), t_f-t_1, N_Simulation, Solver);
X = [X_1, X_2(:,2:end)];
time_vector = [time_vector_1, t_1 + time_vector_2(2:end)];
x_AIF_2 = X(OutputIndex, :);

% AIF Controller 3
Gains.omega_0 = omega_0_vector(3);
Gains.K_I = K_I_vector(3);
Parameters_Disturbed = Parameters;
Parameters_Disturbed.mu = Parameters.mu * Factor_Disturbance;
[Parameters_Disturbed, ~] = Gains2Params(Gains, Parameters_Disturbed, r_Disturbed, SupportingInput);
Parameters = Parameters_Disturbed; Parameters.mu = Parameters_Disturbed.mu / Factor_Disturbance;

IC = FixedPoint(Parameters, SupportingInput);
[time_vector_1, X_1] = DSA(StoichiometryMatrix, PropensityFunction, Parameters, IC, t_1, N_Simulation, Solver);
[time_vector_2, X_2] = DSA(StoichiometryMatrix, PropensityFunction, Parameters_Disturbed, X_1(:,end), t_f-t_1, N_Simulation, Solver);
X = [X_1, X_2(:,2:end)];
time_vector = [time_vector_1, t_1 + time_vector_2(2:end)];
x_AIF_3 = X(OutputIndex, :);
        
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
LineWidth_Thin = ScalingFactor*0.25 * SS;
MarkerSize = ScalingFactor*7 * SS;
PIColor1 = [55,126,184]/255; % Blue
PIColor2 = [77,175,74]/255; % Green
PIColor3 = [228,26,28]/255; % Red


%% Set Figure 1
Figure1_Name = 'AIF_RootLocus_BD_Simulation';
Handle_Figure1 = figure();
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [2, 2, Figure_Width/2, Figure_Height*1.2];
    Handle_Figure1.PaperPositionMode = 'auto';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];

%% Axis for Simulations
Handle_Axis1 = axes(Handle_Figure1);
    Handle_Axis1.Position = [0.15, 0.2, 0.8, 0.7];
    Handle_Axis1.Box = 'on';
    Handle_Axis1.FontSize = FontSize;
    Handle_Axis1.LineWidth = LineWidth_Thin;
    Handle_Axis1.FontSize = FontSize;
    hold(Handle_Axis1, 'on');
    grid(Handle_Axis1, 'on');
    Handle_Axis1.XMinorGrid = 'off';
    Handle_Axis1.YMinorGrid = 'off';
    Handle_Axis1.XLim = [0.97*t_1, t_f];
    Handle_Axis1.XTick = t_1 : 10 : t_f;
    Handle_Axis1.XTickLabel = 0 : 10 : t_f - t_1;
    Handle_Axis1.XLabel.String = 'Time';
    % Handle_Axis1.YLabel.String = 'Output';
    Handle_Axis1.YLim = [0.9*r, 1.1*r_Disturbed];
    Handle_Axis1.Title.String = 'Simulations (Output)';
plot(Handle_Axis1, [0, t_1, t_1, t_f], [r, r, r_Disturbed, r_Disturbed], 'LineWidth', LineWidth_Thin, 'Color', 'k', 'LineStyle', '--');
plot(Handle_Axis1, time_vector, x_AIF_1, 'LineWidth', LineWidth, 'Color', PIColor1);
plot(Handle_Axis1, time_vector, x_AIF_2, 'LineWidth', LineWidth, 'Color', PIColor2);
plot(Handle_Axis1, time_vector, x_AIF_3, 'LineWidth', LineWidth, 'Color', PIColor3);
Handle_Legend = legend(Handle_Axis1, {'Setpoint', ['$K_I = ', num2str(K_I_vector(1)), ', \omega_0 = ', num2str(omega_0_vector(1)), '$'], ['$K_I = ', num2str(K_I_vector(2)), ', \omega_0 = ', num2str(omega_0_vector(2)),  '$'], ['$K_I = ', num2str(K_I_vector(3)), ', \omega_0 = ', num2str(omega_0_vector(3)),  '$']});
Handle_Legend.Location = 'southeast';
Handle_Legend.Interpreter = 'latex';
Handle_Legend.FontSize = FontSize/1.1;

%% Save Figures
if Save_Flag == 1  
    Handle_Figure1.Color = 'none';
    set(Handle_Figure1, 'InvertHardCopy', 'off');
  	print(Handle_Figure1, Figure1_Name, '-dpdf', '-painters','-bestfit');
    Handle_Figure1.Color = [1, 1, 1];
end