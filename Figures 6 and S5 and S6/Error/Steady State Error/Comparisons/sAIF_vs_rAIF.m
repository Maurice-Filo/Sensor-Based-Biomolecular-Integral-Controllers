%% Clear Workspace
close all
clear
clc
disp(['Running: "', mfilename, '.m"']);

%% Fixed Parameters
Parameters.gamma_1 = 0.1;
Parameters.mu = 10;
Parameters.theta = 1;
Parameters.eta = 10;
Parameters.delta = 0.1;
DisturbedParameters = Parameters;

%% Disturbance
Disturbance = 'gamma_1';
Disturbance_Factor = 0.5;
rho = 1 / Parameters.gamma_1;

%% Boundaries
r_Fine = linspace(0, 20, 1000);
f_2 = (1/2/Parameters.eta/Parameters.delta) * (Parameters.eta*(Parameters.theta*r_Fine - Parameters.mu) - Parameters.delta^2 + sqrt( (Parameters.eta*(Parameters.theta*r_Fine - Parameters.mu) - Parameters.delta^2).^2 + 4*Parameters.delta^2*Parameters.eta*Parameters.theta*r_Fine) );
f_1 = Parameters.mu ./ (Parameters.eta * f_2 + Parameters.delta);
G_Bound_rAIF = r_Fine ./ (rho * f_1);
G_Bound_sAIF = r_Fine ./ (rho * f_2);

%% Network Construction
FixedPoint_rAIF = @ComputeFP_rAIF_BirthDeath;
FixedPoint_sAIF = @ComputeFP_sAIF_BirthDeath;
Output_Index = 1;

%% Grid
N_Gains = 300;
N_SetPoints = 320;
Gains_vector = logspace(-6, 5, N_Gains);
SetPoints_vector = linspace(1, 20, N_SetPoints);

%% Simulations for rAIF & sAIF
SimTime = tic;
SS_PreDisturbance_rAIF = zeros(N_SetPoints, N_Gains);
SS_PostDisturbance_rAIF = zeros(N_SetPoints, N_Gains);
SS_PreDisturbance_sAIF = zeros(N_SetPoints, N_Gains);
SS_PostDisturbance_sAIF = zeros(N_SetPoints, N_Gains);
SetPoints_Matrix = zeros(N_SetPoints, N_Gains);
Gains_Matrix = zeros(N_SetPoints, N_Gains);
for i = 1 : N_SetPoints
    r = SetPoints_vector(i);
    disp(['Progress: (i,j) = (', num2str(i), ',:) out of (', num2str(N_SetPoints), ',', num2str(N_Gains), ')']);
    for j = 1 : N_Gains
        G = Gains_vector(j);
        u = r / rho;
        z_2 = (1/2/Parameters.eta/Parameters.delta) * ( Parameters.eta*(Parameters.theta*r - Parameters.mu) - Parameters.delta^2 + sqrt( (Parameters.eta*(Parameters.theta*r - Parameters.mu) - Parameters.delta^2)^2 + 4*Parameters.delta^2*Parameters.eta*Parameters.theta*r) );
        z_1 = Parameters.mu / (Parameters.eta*z_2 + Parameters.delta);
        % rAIF
        Parameters_rAIF = Parameters;
        Parameters_rAIF.alpha = u^2 / (u - G*z_1);
        Parameters_rAIF.kappa = G*z_1^2 / (u - G*z_1);
        DisturbedParameters_rAIF = Parameters_rAIF;
        DisturbedParameters_rAIF.(Disturbance) = Parameters_rAIF.(Disturbance) * Disturbance_Factor;
        FP_rAIF = FixedPoint_rAIF(Parameters_rAIF); 
        SS_PreDisturbance_rAIF(i,j) = FP_rAIF(Output_Index);
        FP_rAIF = FixedPoint_rAIF(DisturbedParameters_rAIF); 
        SS_PostDisturbance_rAIF(i,j) = FP_rAIF(Output_Index);
        % sAIF
        Parameters_sAIF = Parameters;
        Parameters_sAIF.alpha = u^2 / (u - G*z_2);
        Parameters_sAIF.kappa = u/G - z_2;
        DisturbedParameters_sAIF = Parameters_sAIF;
        DisturbedParameters_sAIF.(Disturbance) = Parameters_sAIF.(Disturbance) * Disturbance_Factor;
        FP_sAIF = FixedPoint_sAIF(Parameters_sAIF); 
        SS_PreDisturbance_sAIF(i,j) = FP_sAIF(Output_Index);
        FP_sAIF = FixedPoint_sAIF(DisturbedParameters_sAIF); 
        SS_PostDisturbance_sAIF(i,j) = FP_sAIF(Output_Index);
        SetPoints_Matrix(i,j) = r;
        Gains_Matrix(i,j) = G;
    end
end
SimTime = toc(SimTime);
SS_Error_rAIF = abs(SS_PreDisturbance_rAIF - SS_PostDisturbance_rAIF) ./ SS_PreDisturbance_rAIF;
SS_Error_sAIF = abs(SS_PreDisturbance_sAIF - SS_PostDisturbance_sAIF) ./ SS_PreDisturbance_sAIF;
disp(['Successfully completed the run of: "', mfilename, '.m". Simulation Time is: ', num2str(SimTime/3600), ' hours.']);

%% Figure Settings
SS = 4;
Figure_Width = 7 * SS;
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
                 228,26,28]/255; ... % Red

%% Set Figure 1
Figure1_Name = 'SS_Error_BirthDeath_rAIF_Comparison';
Handle_Figure1 = figure;
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [0, Figure_Height, Figure_Width/2, Figure_Height];
    Handle_Figure1.PaperPositionMode = 'auto';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];

%% Set Axis 1: Topology rAIF
Handle_Axis1 = axes(Handle_Figure1);
    Handle_Axis1.Position = [0.12, 0.14, 0.79, 0.82];
    Handle_Axis1.Box = 'on';
    Handle_Axis1.FontSize = FontSize;
    hold(Handle_Axis1, 'on');
    grid(Handle_Axis1, 'on');
    Handle_Axis1.Layer = 'top';
    Handle_Axis1.XMinorGrid = 'off';
    Handle_Axis1.YMinorGrid = 'off';
    Handle_Axis1.XScale = 'linear';
    Handle_Axis1.YScale = 'log';
    Handle_Axis1.XLabel.Interpreter = 'latex';
    Handle_Axis1.YLabel.Interpreter = 'latex';
    Handle_Axis1.ZScale = 'log';
    % Handle_Axis1.XLim = [5, 15];
    Handle_Axis1.YLim = [1e-3, 10];
    Handle_Axis1.XTick = [0, 5, 10, 15, 20];
    Handle_Axis1.YTick = [1e-3, 1e-2, 1e-1, 1];
    % Handle_Axis1.ZTick = [1e-2, 1e-1, 1];

% Handle_Axis1.XLabel.String = 'Setpoint';
% Handle_Axis1.YLabel.String = 'Gain';
% Handle_Axis1.ZLabel.String = 'SS Error';

plot(Handle_Axis1, r_Fine, G_Bound_rAIF, 'LineWidth', LineWidth, 'Color', NominalColors(1,:));
plot(Handle_Axis1, r_Fine, G_Bound_sAIF, 'LineWidth', LineWidth, 'Color', NominalColors(4,:));
% scatter3(Handle_Axis1, SetPoints_Matrix(:), Gains_Matrix(:), SS_Error_rAIF(:), MarkerSize, [0, 0, 1], 'filled');
% scatter3(Handle_Axis1, SetPoints_Matrix(:), Gains_Matrix(:), SS_Error_sAIF(:), MarkerSize, [1, 0, 0], 'filled');
Handle_rAIF = surf(Handle_Axis1, SetPoints_Matrix, Gains_Matrix, SS_Error_rAIF, 'FaceColor', NominalColors(4,:), 'EdgeColor', 'none');
Handle_sAIF = surf(Handle_Axis1, SetPoints_Matrix, Gains_Matrix, SS_Error_sAIF, 'FaceColor', NominalColors(1,:), 'EdgeColor', 'none');
Handle_rAIF.FaceAlpha = 0.7;
Handle_sAIF.FaceAlpha = 0.7;

N_Skip = 5;
surf(Handle_Axis1, SetPoints_Matrix(1:N_Skip:end, 1:N_Skip:end), Gains_Matrix(1:N_Skip:end, 1:N_Skip:end), SS_Error_rAIF(1:N_Skip:end, 1:N_Skip:end), 'FaceColor', 'none', 'EdgeAlpha', 0.7);
surf(Handle_Axis1, SetPoints_Matrix(1:N_Skip:end, 1:N_Skip:end), Gains_Matrix(1:N_Skip:end, 1:N_Skip:end), SS_Error_sAIF(1:N_Skip:end, 1:N_Skip:end), 'FaceColor', 'none', 'EdgeAlpha', 0.7);

x_coord = Parameters.mu/Parameters.theta - Parameters.delta^2/Parameters.eta/Parameters.theta;
y_range = Handle_Axis1.YLim;
z_range = Handle_Axis1.ZLim;
vertices = [
    x_coord, y_range(1), z_range(1); 
    x_coord, y_range(2), z_range(1);
    x_coord, y_range(2), z_range(2);
    x_coord, y_range(1), z_range(2);
];
faces = [1, 2, 3, 4]; 

patch(Handle_Axis1, 'Vertices', vertices, 'Faces', faces, ...
      'FaceColor', 0.5*ones(1,3), 'FaceAlpha', 0.5, 'EdgeColor', 'none');

% legend(Handle_Axis1, [Handle_rAIF, Handle_sAIF], {'rAIF', 'sAIF'});

%% Revoloving the Plot
Save_Flag = 0;
Handle_Figure1.Color = [1, 1, 1];
% Figure1_Name = 'rAIF_sAIF_Animation.gif';
Figure1_Name = 'rAIF_sAIF';
EZ = 20;
az_vector = linspace(-37.5,-37.5 + 360, 100);
view(Handle_Axis1, -13.0357, 27.3763);
print(Handle_Figure1, Figure1_Name, '-dpng');
% pause()
% for i = 1 : length(az_vector)
%         view(Handle_Axis1, az_vector(i), EZ);
%         drawnow();
%         if Save_Flag == 1
%             frame = getframe(gcf);
%             img = frame2im(frame);
%             [imgIndexed, colormap] = rgb2ind(img, 24); % Convert to indexed image
%             if i == 1
%                 imwrite(imgIndexed, colormap, Figure1_Name, 'gif', 'LoopCount', inf, 'DelayTime', 0.01);
%             else
%                 imwrite(imgIndexed, colormap, Figure1_Name, 'gif', 'WriteMode', 'append', 'DelayTime', 0.01);
%             end
%             % print(Handle_Figure1, [Figure1_Name, num2str(i)], '-dpdf');
%         end
% end
