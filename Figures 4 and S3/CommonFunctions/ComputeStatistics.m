function [t_vector, MeanX, VarX, FinalDistribution] = ComputeStatistics(T, X, N)
% This function takes the stochastic trajectories of multiple species and
% compute the mean trajectory, variance trajectory and final distribution. 
% The time horizon is discretized to N samples, and the stochastic trajectories are zero order
% held between the jumping times.

    %% Find the maximum time among the different stochastic trajectories
    t_max = 1e20;
    for i = 1 : length(T)
        T_Vector = T{i};
        t_max = min(t_max, T_Vector(end));
    end
    
    %% Discretize the Time Horizon
    t_vector = linspace(0, t_max, N);
        
    %% Compute the sample paths at the common time horizon
    X_Resampled = cell(length(T), 1);
    for i = 1 : length(T)
        [X_Resampled{i}] = ResampleZOH(X{i}, T{i}, t_vector);
    end
    %% Compute the Statistics
    MeanX = zeros(size(X{1}, 1), length(t_vector));
    VarX = zeros(size(X{1}, 1), length(t_vector));
    for j = 1 : size(X{1}, 1)
        Xj = zeros(length(T), length(t_vector));
        for i = 1 : length(T)
            Xj(i,:) = X_Resampled{i}(j,:);
        end
        MeanX(j,:) = mean(Xj, 1);
        VarX(j,:) = var(Xj, 1);
    end
    MeanX(:,end) = [];
    VarX(:,end) = [];
    t_vector(end) = [];
    
    %% Compute Final Distribution
    FinalDistribution = zeros(size(X{1},1), length(T));
  	for i = 1 : length(T)
        FinalDistribution(:,i) = X{i}(:,end);
  	end
end

