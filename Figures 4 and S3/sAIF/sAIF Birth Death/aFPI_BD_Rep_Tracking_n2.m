%% Clear Workspace
close all
clear
clc

Save_Flag = 0;

%% Functions
StoichiometryMatrix = StoichiometryMatrix_aFPI_Rep_BD();
PropensityFunction = @PropensityFunction_aFPI_Rep_BD;
FixedPoint = @aFPI_Rep_FixedPoint;
Params2Gains = @aFPI_Rep_Params2Gains;
Gains2Params = @aFPI_Rep_Gains2Params;
SupportingInput = @SupportingInput_BD;
OutputIndex = 1;

%% Plant Parameters
Parameters.gamma_1 = 1;

%% Controller Parameters
Parameters.n = 2;
Parameters.mu = 5;
r = 5;
Parameters.theta = Parameters.mu / r;
K_S = Parameters.theta;

%% Disturbance
Tracking_Factor = 2;
r_New = r * Tracking_Factor;
Parameters_New = Parameters;
Parameters_New.mu = Parameters.mu * Tracking_Factor;
u_bar = SupportingInput(Parameters, r_New);

%% Simulation Settings
N_Simulation = 1000;
Solver = 'ODE15s';
t_1 = 100;
t_f = 120;

%% Sweep the Pole
N = 1000;
N_Nominal = 4;
a_Sweep = linspace(0.8*Parameters.gamma_1, 2.1*Parameters.gamma_1, N+1);
a_Sweep = a_Sweep(2:end); % to avoid numerical error when a very close to gamma_1/2
Indeces_Nominal = [1, 200, 500, 1000];

%% Compute Gains and Biological Parameters for each placed pole
Gains = cell(N,1);
omega_c_vector = zeros(N,1);
K_P_vector = zeros(N,1);
K_I_vector = zeros(N,1);
CoverageFlag = zeros(N,1);
BioParameters = cell(N,1);
Feasibility_Flag = zeros(N,1);
eta_vector = zeros(N,1);
kappa_vector = zeros(N,1);
alpha_vector = zeros(N,1);
for i = 1 : N
    a = a_Sweep(i);
    omega_c = 3*a - Parameters.gamma_1;
    K_I = a^3 / (K_S * omega_c);
    K_P = (3*a^2-Parameters.gamma_1*omega_c)/(K_S*omega_c);
    Gains{i}.K_P = K_P;
    Gains{i}.K_I = K_I;
    Gains{i}.omega_c = omega_c;
    Gains{i}.K_S = K_S; 
    omega_c_vector(i) = omega_c;
    K_I_vector(i) = K_I;
    K_P_vector(i) = K_P;
    CoverageFlag(i) = (K_P >= 0) && (K_I >= 0) && (omega_c >= 0);
    [BioParameters{i}, Feasibility_Flag(i)] = Gains2Params(Gains{i}, Parameters_New, r_New, SupportingInput);
    kappa_vector(i) = BioParameters{i}.kappa;
    eta_vector(i) = BioParameters{i}.eta;
    alpha_vector(i) = BioParameters{i}.alpha;
end

%% Simulations
a_Nominal = zeros(N_Nominal,1);
Gains_Nominal = cell(N_Nominal,1);
BioParameters_Nominal = cell(N_Nominal,1);
x_FPI = zeros(N_Nominal, 2*N_Simulation-1);
for i = 1 : N_Nominal
    Index = Indeces_Nominal(i);
    a_Nominal(i) = a_Sweep(Index);
    Parameters = BioParameters{Index};
    Parameters.mu = Parameters.mu / Tracking_Factor;
    IC = FixedPoint(Parameters, SupportingInput);
    [time_vector_1, X_1] = DSA(StoichiometryMatrix, PropensityFunction, Parameters, IC, t_1, N_Simulation, Solver);
    [time_vector_2, X_2] = DSA(StoichiometryMatrix, PropensityFunction, BioParameters{Index}, X_1(:,end), t_f-t_1, N_Simulation, Solver);
    X = [X_1, X_2(:,2:end)];
    time_vector = [time_vector_1, t_1 + time_vector_2(2:end)];
    x_FPI(i,:) = X(OutputIndex,:);
    Gains_Nominal{i} = Params2Gains(BioParameters{Index}, SupportingInput);
    BioParameters_Nominal{i} = BioParameters{Index};
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
NominalColors = [55,126,184; ... % Blue
                 77,175,74; ... % Green
                 255,127,0; ... % Orange
                 228,26,28]/255; ... % Red
             
             
%% Set Figure 1
Figure1_Name = 'FPI_Rep_BD_Poles_n2';
Handle_Figure1 = figure();
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [2, 2, Figure_Width/2, Figure_Height*1.2];
    Handle_Figure1.PaperPositionMode = 'auto';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];
    
%% Set Figure 2
Figure2_Name = 'FPI_Rep_BD_Gains_n2';
Handle_Figure2 = figure();
    Handle_Figure2.PaperUnits = 'centimeters';
    Handle_Figure2.Units = 'centimeters';
    Handle_Figure2.Position = [2, 2, Figure_Width/2, Figure_Height*1.2];
    Handle_Figure2.PaperPositionMode = 'auto';
    Handle_Figure2.PaperSize = [Handle_Figure2.PaperPosition(3), Handle_Figure2.PaperPosition(4)];
    
%% Set Figure 3
Figure3_Name = 'FPI_Rep_BD_BioParameters_n2';
Handle_Figure3 = figure();
    Handle_Figure3.PaperUnits = 'centimeters';
    Handle_Figure3.Units = 'centimeters';
    Handle_Figure3.Position = [2, 2, Figure_Width/2, Figure_Height*1.2];
    Handle_Figure3.PaperPositionMode = 'auto';
    Handle_Figure3.PaperSize = [Handle_Figure3.PaperPosition(3), Handle_Figure3.PaperPosition(4)];

%% Set Figure 4
Figure4_Name = 'FPI_Rep_BD_Simulations_n2';
Handle_Figure4 = figure();
    Handle_Figure4.PaperUnits = 'centimeters';
    Handle_Figure4.Units = 'centimeters';
    Handle_Figure4.Position = [2, 2, Figure_Width/2, Figure_Height*1.2];
    Handle_Figure4.PaperPositionMode = 'auto';
    Handle_Figure4.PaperSize = [Handle_Figure4.PaperPosition(3), Handle_Figure4.PaperPosition(4)];

%% Axis for Poles
Handle_Axis1 = axes(Handle_Figure1);
    Handle_Axis1.Position = [0.05, 0.16, 0.78, 0.76];
    Handle_Axis1.Box = 'on';
    Handle_Axis1.BoxStyle = 'full';
    Handle_Axis1.LineWidth = LineWidth_Thin;
    Handle_Axis1.FontSize = FontSize_Small*1.25;
    hold(Handle_Axis1, 'on');
    grid(Handle_Axis1, 'on');
    Handle_Axis1.XMinorGrid = 'off';
    Handle_Axis1.YMinorGrid = 'off';
    Handle_Axis1.XLim = [-10, 0.5];
    Handle_Axis1.XTick = [-10, -8, -6, -4, -2, -Parameters.gamma_1, 0];
    Handle_Axis1.XTickLabel = {-10, -8, -6, -4, -2, '$-\gamma_1$', 0};
    Handle_Axis1.TickLabelInterpreter = 'latex';
%     Handle_Axis1.XTick = [-Parameters.gamma_1*2.1769, -Parameters.gamma_1, -Parameters.gamma_1*0.7082, 0];
    Handle_Axis1.YLim = [-0.4, 0.4];
    Handle_Axis1.YTick = [-0.2, 0, 0.2];
    Handle_Axis1.XLabel.String = 'Real Axis';
    Handle_Axis1.XLabel.FontSize = FontSize;
    Handle_Axis1.YLabel.String = 'Imaginary Axis';
    Handle_Axis1.YLabel.FontSize = FontSize;
    colormap(Handle_Axis1, 'turbo');
    Handle_Axis1.YAxisLocation = 'Right';
%     Handle_Axis1.Title.String = 'Pole Placement';
    Handle_Axis1.Title.FontSize = FontSize;
scatter(Handle_Axis1, -a_Sweep, 0*a_Sweep, MarkerSize*3, a_Sweep, 'filled');
for i = 1 : N_Nominal
    plot(Handle_Axis1, -a_Nominal(i), 0, 'Marker', 'o', 'MarkerSize', MarkerSize_Small, 'Color', 'k', 'LineStyle', 'none', 'MarkerFaceColor', NominalColors(i,:));
end
% plot(Handle_Axis1, [0, 0], Handle_Axis1.YLim, 'LineWidth', LineWidth_Thin, 'Color', 'k', 'LineStyle', ':');
% xline(Handle_Axis1, -Parameters.gamma_1*0.7082, 'LineWidth', LineWidth_Thin, 'Color', 'k', 'LineStyle', '--');
% xline(Handle_Axis1, -Parameters.gamma_1*2.1769, 'LineWidth', LineWidth_Thin, 'Color', 'k', 'LineStyle', '--');
    
%% Axis for Gains
Handle_Axis2 = axes(Handle_Figure2);
    Handle_Axis2.Position = [0.18, 0.14, 0.7, 0.78];
    Handle_Axis2.Box = 'on';
    Handle_Axis2.BoxStyle = 'full';
    Handle_Axis2.LineWidth = LineWidth_Thin;
    Handle_Axis2.FontSize = FontSize;
    hold(Handle_Axis2, 'on');
    grid(Handle_Axis2, 'on');
    Handle_Axis2.XMinorGrid = 'off';
    Handle_Axis2.YMinorGrid = 'off';
    Handle_Axis2.XLabel.String = '$K_P$';
    Handle_Axis2.YLabel.String = '$K_I$';
    Handle_Axis2.ZLabel.String = '$\omega_0$';
    Handle_Axis2.XLabel.Interpreter = 'latex';
    Handle_Axis2.YLabel.Interpreter = 'latex';
    Handle_Axis2.ZLabel.Interpreter = 'latex';
    view(Handle_Axis2, [40, 25])
    colormap(Handle_Axis2, 'turbo');
%  	Handle_Axis2.Title.String = 'PI Gains & Cutoff Frequency';
    Handle_Axis2.Title.FontSize = FontSize*1;
    Handle_Axis2.XTick = [0.4, 0.8, 1.2];
    Handle_Axis2.YTick = [0.5, 1, 1.5];
scatter3(Handle_Axis2, K_P_vector, K_I_vector, omega_c_vector, MarkerSize*3, linspace(0, 1, N), 'filled');
for i = 1 : N_Nominal
    plot3(Handle_Axis2, Gains_Nominal{i}.K_P, Gains_Nominal{i}.K_I, Gains_Nominal{i}.omega_c, 'Marker', 'o', 'MarkerSize', MarkerSize, 'Color', 'k', 'LineStyle', 'none', 'MarkerFaceColor', NominalColors(i,:));
end

%% Axis for BioParameters
Handle_Axis3 = axes(Handle_Figure3);
    Handle_Axis3.Position = [0.18, 0.13, 0.72, 0.75];
    Handle_Axis3.Box = 'on';
    Handle_Axis3.BoxStyle = 'full';
    Handle_Axis3.LineWidth = LineWidth_Thin;
    Handle_Axis3.FontSize = FontSize;
    hold(Handle_Axis3, 'on');
    grid(Handle_Axis3, 'on');
    Handle_Axis3.XMinorGrid = 'off';
    Handle_Axis3.YMinorGrid = 'off';
    Handle_Axis3.XLabel.String = '$\gamma$';
    Handle_Axis3.YLabel.String = '$\alpha$';
    Handle_Axis3.ZLabel.String = '$\eta$';
    Handle_Axis3.XLabel.Interpreter = 'latex';
    Handle_Axis3.YLabel.Interpreter = 'latex';
    Handle_Axis3.ZLabel.Interpreter = 'latex';
    view(Handle_Axis3, [40, 30])
    colormap(Handle_Axis3, 'turbo');
%     Handle_Axis3.Title.String = 'Biomolecular Parameters';
    Handle_Axis3.Title.FontSize = FontSize*0.95;
scatter3(Handle_Axis3, kappa_vector, alpha_vector, eta_vector, MarkerSize*3, linspace(0, 1, N), 'filled');
for i = 1 : N_Nominal
    plot3(Handle_Axis3, BioParameters_Nominal{i}.kappa, BioParameters_Nominal{i}.alpha, BioParameters_Nominal{i}.eta, 'Marker', 'o', 'MarkerSize', MarkerSize, 'Color', 'k', 'LineStyle', 'none', 'MarkerFaceColor', NominalColors(i,:));
end

%% Axis for Simulations
Handle_Axis4 = axes(Handle_Figure4);
    Handle_Axis4.Position = [0.15, 0.2, 0.8, 0.7];
    Handle_Axis4.Box = 'on';
    Handle_Axis4.FontSize = FontSize;
    Handle_Axis4.LineWidth = LineWidth_Thin;
    Handle_Axis4.FontSize = FontSize;
    hold(Handle_Axis4, 'on');
    grid(Handle_Axis4, 'on');
    Handle_Axis4.XMinorGrid = 'off';
    Handle_Axis4.YMinorGrid = 'off';
    Handle_Axis4.XLim = [0.97*t_1, t_f];
    Handle_Axis4.XTick = t_1 : 10 : t_f;
    Handle_Axis4.XTickLabel = 0 : 10 : t_f - t_1;
    Handle_Axis4.XLabel.String = 'Time';
    % Handle_Axis4.YLabel.String = 'Output';
    Handle_Axis4.YLim = [0.9*r, 1.1*r_New];
    % Handle_Axis4.Title.String = 'Simulations';
plot(Handle_Axis4, [0, t_1, t_1, t_f], [r, r, r_New, r_New], 'LineWidth', LineWidth_Thin, 'Color', 'k', 'LineStyle', '--');
for i = 1 : N_Nominal
    plot(Handle_Axis4, time_vector, x_FPI(i,:), 'LineWidth', LineWidth, 'Color', NominalColors(i,:));
end

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
    Handle_Figure4.Color = 'none';
    set(Handle_Figure4, 'InvertHardCopy', 'off');
    print(Handle_Figure4, Figure4_Name, '-dpdf', '-painters','-bestfit');
    Handle_Figure4.Color = [1, 1, 1];
end

