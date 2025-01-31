import pandas as pd
import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
import pdb

def Visualize_Dynamics(t, x, y, Color_x, Color_y, Opacity, Title, Settings, Axis=None):
    # Compute Means
    Mean_x = x.mean(1)
    Mean_y = y.mean(1)
    spline_x = UnivariateSpline(t, Mean_x, s=0) 
    spline_y = UnivariateSpline(t, Mean_y, s=0)
    Mean_x = spline_x(t)
    Mean_y = spline_y(t)
    # Compute Standard Deviations
    Std_x = x.std(1)
    Std_y = y.std(1)
    # Create Figure and Axis
    if Axis is None:
        Figure, Axis = plt.subplots(figsize=Settings.get('FigureSize', (6.4, 4)))
    # Disturbance
    Axis.axvline(x=1, color=Settings.get('Color_Disturbance', 'k'), linestyle='--')
    # Data Points
    Axis.plot(t, x, marker='o', linestyle='', markeredgecolor='none', markerfacecolor=Color_x, markersize=Settings.get('MarkerSize', 3))
    Axis.plot(t, y, marker='o', linestyle='', markeredgecolor='none', markerfacecolor=Color_y, markersize=Settings.get('MarkerSize', 3))
    # Plot Means
    Axis.plot(t, Mean_x, color=Color_x, label='Without Disturbance')
    Axis.plot(t, Mean_y, color=Color_y, label='With Disturbance')
    # Plot Standard Deviations
    Axis.fill_between(t, Mean_x-Std_x, Mean_x+Std_x, color=Color_x, alpha=Opacity, edgecolor='none')
    Axis.fill_between(t, Mean_y-Std_y, Mean_y+Std_y, color=Color_y, alpha=Opacity, edgecolor='none')
    # Axis Limits
    if 'XLimits' in Settings:
        Axis.set_xlim(Settings['XLimits'])
    if 'YLimits' in Settings:
        Axis.set_ylim(Settings['YLimits'])
    # Title
    Axis.set_title(Title, fontweight='bold')
    # Axis Labels
    Axis.set_xlabel('Time (h)', fontsize = Settings.get('FontSize', 12))
    Axis.set_ylabel('Mean mScarlet-I-H (MEF)', fontsize = Settings.get('FontSize', 12))
    # Ticks Fontsize
    Axis.tick_params(axis='both', which='major', labelsize=Settings.get('FontSize', 12))
    # Grid
    Axis.set_axisbelow(True)
    Axis.grid(True)
    Axis.grid(color='gray', linestyle='--', linewidth=0.5)
    # Fit Axis to Figure
    plt.tight_layout()    
    # Return Handles
    Figure = Axis.figure
    return Figure, Axis

def Visualize_CV_Mean(mu, cv, Color, Opacity, Title, Settings, Process, Axis=None):
    # Create Figure and Axis
    if Axis is None:
        Figure, Axis = plt.subplots(figsize=Settings.get('FigureSize', (6.4, 4)))
    # Data Points
    Points = Axis.plot(mu, cv, marker='o', linestyle='', markeredgecolor='black', markeredgewidth=0.01, markerfacecolor=Color, markersize=Settings.get('MarkerSize', 2), alpha=Opacity)
    if Process:
        # Compute Means
        Mean_mu = mu.mean(1)
        Mean_cv = cv.mean(1)
        # Compute Standard Deviations
        Std_mu = mu.std(1)
        Std_cv = cv.std(1)
        # Plot Means
        Axis.plot(Mean_mu, Mean_cv, color=Color, linestyle='--', alpha=Opacity)
        # Error Bars
        Axis.errorbar(Mean_mu, Mean_cv, xerr=Std_mu, yerr=Std_cv, fmt='o', color=Color, capsize=4, ecolor=Color, linestyle='None', markersize=2*Settings.get('MarkerSize', 2), alpha=Opacity)
    # Axis Limits
    if 'XLimits' in Settings:
        Axis.set_xlim(Settings['XLimits'])
    if 'YLimits' in Settings:
        Axis.set_ylim(Settings['YLimits'])
    # Title
    Axis.set_title(Title, fontweight='bold')
    # Ticks Fontsize
    Axis.tick_params(axis='both', which='major', labelsize=Settings.get('FontSize', 12))
    # Axis Labels
    Axis.set_xlabel('Mean (MEF)', fontsize = Settings.get('FontSize', 12))
    Axis.set_ylabel('CV', fontsize = Settings.get('FontSize', 12))
    # Grid
    Axis.set_axisbelow(True)
    Axis.grid(True)
    Axis.grid(color='gray', linestyle='--', linewidth=0.5)
    # Fit Axis to Figure
    plt.tight_layout() 
    # Return Handles
    Figure = Axis.figure
    return Figure, Axis, Points

def Visualize_RPA1(Input, Output, Output_Disturbance, Color, Color_Disturbance, Width, Title, Settings, Normalize, XLabels, Axis=None):
    # Create Figure and Axis
    if Axis is None:
        Figure, Axis = plt.subplots(figsize=Settings.get('FigureSize', (8, 4)))
    # Compute Means
    Mean = Output.mean(1)
    Mean_Disturbance = Output_Disturbance.mean(1)
    # Compute Standard Deviations
    Std = Output.std(1)
    Std_Disturbance = Output_Disturbance.std(1)
    # Normalize
    if Normalize:
        Std_Disturbance = Std_Disturbance / Mean
        Std = Std / Mean
        Mean_Disturbance = Mean_Disturbance / Mean
        Mean = Mean / Mean
    # Plot Means
    Axis.bar(np.arange(1, Input.size+1) - Width/2 - 0.3*Width, Mean, Width, color=Color, edgecolor='black')
    Axis.bar(np.arange(1, Input.size+1) + Width/2 + 0.3*Width, Mean_Disturbance, Width, color=Color_Disturbance, edgecolor='black')
    # Error Bars
    Axis.errorbar(np.arange(1, Input.size+1) - Width/2 - 0.3*Width, Mean, yerr=Std, fmt='None', capsize=4, ecolor='black', linestyle='None')
    Axis.errorbar(np.arange(1, Input.size+1) + Width/2 + 0.3*Width, Mean_Disturbance, yerr=Std_Disturbance, fmt='None', capsize=4, ecolor='black', linestyle='None')
    # Plot Data Points
    if Normalize:
        Axis.plot(np.arange(1, Input.size+1) - Width/2 - 0.3*Width, Output.div(Output.mean(1), axis=0), marker='o', linestyle='', markeredgecolor='black', markerfacecolor=Color, markersize=Settings.get('MarkerSize', 2), markeredgewidth=0.2)
        Axis.plot(np.arange(1, Input.size+1) + Width/2 + 0.3*Width, Output_Disturbance.div(Output.mean(1), axis=0), marker='o', linestyle='', markeredgecolor='black', markerfacecolor=Color_Disturbance, markersize=Settings.get('MarkerSize', 2), markeredgewidth=0.2)
    else:
        Axis.plot(np.arange(1, Input.size+1) - Width/2 - 0.3*Width, Output, marker='o', linestyle='', markeredgecolor='black', markerfacecolor=Color, markersize=Settings.get('MarkerSize', 2), markeredgewidth=0.2)        
        Axis.plot(np.arange(1, Input.size+1) + Width/2 + 0.3*Width, Output_Disturbance, marker='o', linestyle='', markeredgecolor='black', markerfacecolor=Color_Disturbance, markersize=Settings.get('MarkerSize', 2), markeredgewidth=0.2)                
    # Axis Limits
    if 'XLimits' in Settings:
        Axis.set_xlim(Settings['XLimits'])
    if 'YLimits' in Settings:
        Axis.set_ylim(Settings['YLimits'])
    # Title
    Axis.set_title(Title, fontweight='bold')
    # Axis Labels
    if XLabels:
        Axis.set_xlabel('Arabinose (%)', fontsize = Settings.get('FontSize', 12))
    if Normalize:
        Axis.set_ylabel('Normalized Output', fontsize = Settings.get('FontSize', 12))
    else:
        Axis.set_ylabel('Mean Output (MEF)', fontsize = Settings.get('FontSize', 12))
    # Axis Ticks
    Axis.set_xticks(np.arange(1, Input.size+1))
    # Axis.set_xticklabels(Input)
    if not XLabels:
        Axis.set_xticklabels([])
    # Ticks Fontsize
    Axis.tick_params(axis='both', which='major', labelsize=Settings.get('FontSize', 12))
    # Grid
    Axis.set_axisbelow(True)
    Axis.grid(True)
    Axis.grid(color='gray', linestyle='--', linewidth=0.5)
    # Fit Axis to Figure
    plt.tight_layout()   
    # Return Handles
    Figure = Axis.figure
    return Figure, Axis

def Visualize_RPA2(Input, Output_1, Output_1_Disturbance, Output_2, Output_2_Disturbance, Color, Color_Disturbance, Width, Title, Legend, Settings, Normalize, XLabels, Axis=None):
    # Create Figure and Axis
    if Axis is None:
        Figure, Axis = plt.subplots(figsize=Settings.get('FigureSize', (8, 4)))
    # Compute Means
    Mean_1 = Output_1.mean(1)
    Mean_1_Disturbance = Output_1_Disturbance.mean(1)
    Mean_2 = Output_2.mean(1)
    Mean_2_Disturbance = Output_2_Disturbance.mean(1)
    # Compute Standard Deviations
    Std_1 = Output_1.std(1)
    Std_1_Disturbance = Output_1_Disturbance.std(1)
    Std_2 = Output_2.std(1)
    Std_2_Disturbance = Output_2_Disturbance.std(1)
    # Normalize
    if Normalize:
        Std_1_Disturbance = Std_1_Disturbance / Mean_1
        Std_1 = Std_1 / Mean_1
        Std_2_Disturbance = Std_2_Disturbance / Mean_2
        Std_2 = Std_2 / Mean_2
        Mean_1_Disturbance = Mean_1_Disturbance / Mean_1
        Mean_1 = Mean_1 / Mean_1
        Mean_2_Disturbance = Mean_2_Disturbance / Mean_2
        Mean_2 = Mean_2 / Mean_2
    # Plot Means
    if Normalize:
        Handle_Mean1 = Axis.bar(np.arange(1, Input.size+1) - Width - 0.3*Width, Mean_1, Width, color=Color[0], edgecolor='black')
        Handle_Mean1_Disturbance = Axis.bar(np.arange(1, Input.size+1), Mean_1_Disturbance, Width, color=Color_Disturbance[0], edgecolor='black')
        Handle_Mean2_Disturbance = Axis.bar(np.arange(1, Input.size+1) + Width + 0.3*Width, Mean_2_Disturbance, Width, color=Color_Disturbance[1], edgecolor='black')
    else:
        Handle_Mean1 = Axis.bar(np.arange(1, Input.size+1) - Width*3/2 - 0.4*Width, Mean_1, Width, color=Color[0], edgecolor='black')
        Handle_Mean1_Disturbance = Axis.bar(np.arange(1, Input.size+1) - Width/2 - 0.2*Width, Mean_1_Disturbance, Width, color=Color_Disturbance[0], edgecolor='black')
        Handle_Mean2 = Axis.bar(np.arange(1, Input.size+1) + Width/2 + 0.2*Width, Mean_2, Width, color=Color[1], edgecolor='black')
        Handle_Mean2_Disturbance = Axis.bar(np.arange(1, Input.size+1) + Width*3/2 + 0.4*Width, Mean_2_Disturbance, Width, color=Color_Disturbance[1], edgecolor='black')
    # Error Bars
    if Normalize:
        Axis.errorbar(np.arange(1, Input.size+1) - Width - 0.3*Width, Mean_1, yerr=Std_1, fmt='None', capsize=4, ecolor='black', linestyle='None')
        Axis.errorbar(np.arange(1, Input.size+1), Mean_1_Disturbance, yerr=Std_1_Disturbance, fmt='None', capsize=4, ecolor='black', linestyle='None')
        Axis.errorbar(np.arange(1, Input.size+1) + Width + 0.3*Width, Mean_2_Disturbance, yerr=Std_2_Disturbance, fmt='None', capsize=4, ecolor='black', linestyle='None')
    else:
        Axis.errorbar(np.arange(1, Input.size+1) - Width*3/2 - 0.4*Width, Mean_1, yerr=Std_1, fmt='None', capsize=4, ecolor='black', linestyle='None')
        Axis.errorbar(np.arange(1, Input.size+1) - Width/2 - 0.2*Width, Mean_1_Disturbance, yerr=Std_1_Disturbance, fmt='None', capsize=4, ecolor='black', linestyle='None')
        Axis.errorbar(np.arange(1, Input.size+1) + Width/2 + 0.2*Width, Mean_2, yerr=Std_2, fmt='None', capsize=4, ecolor='black', linestyle='None')
        Axis.errorbar(np.arange(1, Input.size+1) + Width*3/2 + 0.4*Width, Mean_2_Disturbance, yerr=Std_2_Disturbance, fmt='None', capsize=4, ecolor='black', linestyle='None')
    # Plot Data Points
    if Normalize:
        Axis.plot(np.arange(1, Input.size+1) - Width - 0.3*Width, Output_1.div(Output_1.mean(1), axis=0), marker='o', linestyle='', markeredgecolor='black', markerfacecolor=Color[0], markersize=Settings.get('MarkerSize', 2), markeredgewidth=0.2)
        Axis.plot(np.arange(1, Input.size+1), Output_1_Disturbance.div(Output_1.mean(1), axis=0), marker='o', linestyle='', markeredgecolor='black', markerfacecolor=Color_Disturbance[0], markersize=Settings.get('MarkerSize', 2), markeredgewidth=0.2)
        Axis.plot(np.arange(1, Input.size+1) + Width + 0.3*Width, Output_2_Disturbance.div(Output_2.mean(1), axis=0), marker='o', linestyle='', markeredgecolor='black', markerfacecolor=Color_Disturbance[1], markersize=Settings.get('MarkerSize', 2), markeredgewidth=0.2)
    else:
        Axis.plot(np.arange(1, Input.size+1) - Width*3/2 - 0.4*Width, Output_1, marker='o', linestyle='', markeredgecolor='black', markerfacecolor=Color[0], markersize=Settings.get('MarkerSize', 2), markeredgewidth=0.2)        
        Axis.plot(np.arange(1, Input.size+1) - Width/2 - 0.2*Width, Output_1_Disturbance, marker='o', linestyle='', markeredgecolor='black', markerfacecolor=Color_Disturbance[0], markersize=Settings.get('MarkerSize', 2), markeredgewidth=0.2)        
        Axis.plot(np.arange(1, Input.size+1) + Width/2 + 0.2*Width, Output_2, marker='o', linestyle='', markeredgecolor='black', markerfacecolor=Color[1], markersize=Settings.get('MarkerSize', 2), markeredgewidth=0.2)        
        Axis.plot(np.arange(1, Input.size+1) + Width*3/2 + 0.4*Width, Output_2_Disturbance, marker='o', linestyle='', markeredgecolor='black', markerfacecolor=Color_Disturbance[1], markersize=Settings.get('MarkerSize', 2), markeredgewidth=0.2)          
    # Axis Limits
    if 'XLimits' in Settings:
        Axis.set_xlim(Settings['XLimits'])
    if 'YLimits' in Settings:
        Axis.set_ylim(Settings['YLimits'])
    # Title
    Axis.set_title(Title, fontweight='bold')
    # Axis Labels
    Axis.set_xlabel('Arabinose (%)', fontsize = Settings.get('FontSize', 12))
    if Normalize:
        Axis.set_ylabel('Normalized Output', fontsize = Settings.get('FontSize', 12))
    else:
        Axis.set_ylabel('Mean Output (MEF)', fontsize = Settings.get('FontSize', 12))
    # Axis Ticks
    Axis.set_xticks(np.arange(1, Input.size+1))
    Axis.set_xticklabels(Input)
    if not XLabels:
        Axis.set_xticklabels([])
    # Ticks Fontsize
    Axis.tick_params(axis='both', which='major', labelsize=Settings.get('FontSize', 12))
    # Grid
    Axis.set_axisbelow(True)
    Axis.grid(True)
    Axis.grid(color='gray', linestyle='--', linewidth=0.5)
    # Legend
    if Normalize:
        Axis.legend([Handle_Mean1, Handle_Mean1_Disturbance, Handle_Mean2_Disturbance], ['Without Disturbance', Legend[1], Legend[2]])
    else:
        Axis.legend([Handle_Mean1, Handle_Mean1_Disturbance, Handle_Mean2, Handle_Mean2_Disturbance], Legend)
    # Fit Axis to Figure
    plt.tight_layout()   
    # Return Handles
    Figure = Axis.figure
    return Figure, Axis

def Visualize_RPA1_Sorted(Threshold, Type, Output, Output_Disturbance, Color, Color_Disturbance, Width, Title, Legend, Settings, Normalize, XLabels, Axis=None):
    # Create Figure and Axis
    if Axis is None:
        Figure, Axis = plt.subplots(figsize=Settings.get('FigureSize', (8, 4)))
    # Compute Means
    Mean = Output.mean(1)
    Mean_Disturbance = Output_Disturbance.mean(1)
    # Compute Standard Deviations
    Std = Output.std(1)
    Std_Disturbance = Output_Disturbance.std(1)
    # Filter
    Filtering_Flags = Mean <= Threshold
    Filtering_Flags = Filtering_Flags.values
    Type = Type.loc[Filtering_Flags]
    Output = Output.loc[Filtering_Flags]
    Output_Disturbance = Output_Disturbance.loc[Filtering_Flags]
    Mean = Mean.loc[Filtering_Flags]
    Mean_Disturbance = Mean_Disturbance.loc[Filtering_Flags]
    Std = Std.loc[Filtering_Flags]
    Std_Disturbance = Std_Disturbance.loc[Filtering_Flags]
    # Sort
    Sorting_Indeces = Mean.argsort()
    Type = Type.iloc[Sorting_Indeces]
    Output = Output.iloc[Sorting_Indeces]
    Output_Disturbance = Output_Disturbance.iloc[Sorting_Indeces]
    Mean = Mean.iloc[Sorting_Indeces]
    Mean_Disturbance = Mean_Disturbance.iloc[Sorting_Indeces]
    Std = Std.iloc[Sorting_Indeces]
    Std_Disturbance = Std_Disturbance.iloc[Sorting_Indeces]
    # Print Range
    print(f"Output Range: [{Mean.iloc[0]}, {Mean.iloc[-1]}],  Fold Change: {Mean.iloc[-1]/Mean.iloc[0]}.")
    # Normalize
    Normalize
    if Normalize:
        Output = Output.div(Mean, axis=0)
        Output_Disturbance = Output_Disturbance.div(Mean, axis=0)
        Output.div(Mean, axis=0)
        Std_Disturbance = Std_Disturbance / Mean
        Std = Std / Mean
        Mean_Disturbance = Mean_Disturbance / Mean
        Mean = Mean / Mean
    # Plot Means
    Sorted_Colors = np.where(Type == 1, Color[0], Color[1]) 
    Sorted_Colors_Disturbance = np.where(Type == 1, Color_Disturbance[0], Color_Disturbance[1]) 
    Handle_Mean1 = Axis.bar(np.arange(1, Type.size+1) - Width/2 - Width/10, Mean, Width, color=Sorted_Colors, edgecolor='black')
    Handle_Mean1_Disturbance = Axis.bar(np.arange(1, Type.size+1) + Width/2 + Width/10, Mean_Disturbance, Width, color=Sorted_Colors_Disturbance, edgecolor='black')
    # Error Bars
    Axis.errorbar(np.arange(1, Type.size+1) - Width/2 - Width/10, Mean, yerr=Std, fmt='None', capsize=4, ecolor='black', linestyle='None')
    Axis.errorbar(np.arange(1, Type.size+1) + Width/2 + Width/10, Mean_Disturbance, yerr=Std_Disturbance, fmt='None', capsize=4, ecolor='black', linestyle='None')
    # Plot Data Points
    Axis.plot(np.arange(1, Type.size+1) - Width/2 - Width/10, Output, marker='o', linestyle='', markeredgecolor='black', markerfacecolor=Color[0], markersize=Settings.get('MarkerSize', 2), markeredgewidth=0.2)        
    Axis.plot(np.arange(1, Type.size+1) + Width/2 + Width/10, Output_Disturbance, marker='o', linestyle='', markeredgecolor='black', markerfacecolor=Color_Disturbance[0], markersize=Settings.get('MarkerSize', 2), markeredgewidth=0.2)        
    # Axis Limits
    if 'XLimits' in Settings:
        Axis.set_xlim(Settings['XLimits'])
    if 'YLimits' in Settings:
        Axis.set_ylim(Settings['YLimits'])
    # Title
    Axis.set_title(Title, fontweight='bold')
    # Axis Labels
    if Normalize:
        Axis.set_ylabel('Normalized Output', fontsize = Settings.get('FontSize', 12))
    else:
        Axis.set_ylabel('Mean Output (MEF)', fontsize = Settings.get('FontSize', 12))
    # Axis Ticks
    Axis.set_xticks(np.arange(1, Type.size+1))
    # Axis.set_xticklabels(Type)
    if not XLabels:
        Axis.set_xticklabels([])
    # Ticks Fontsize
    Axis.tick_params(axis='both', which='major', labelsize=Settings.get('FontSize', 12))
    # Grid
    Axis.set_axisbelow(True)
    Axis.grid(True)
    Axis.grid(color='gray', linestyle='--', linewidth=0.5)
    # Average Error
    if Normalize: 
        for i, value in enumerate(Mean_Disturbance):
            Axis.text(i+1, Mean_Disturbance.iloc[i] * 1.1, str(round((Mean_Disturbance.iloc[i]-1)*100))+ '% Error', ha='center', va='bottom', color='r', fontsize=Settings.get('FontSize', 12))
    # if Normalize:
    #     Avg_Disturbance = Mean_Disturbance.mean()
    #     Axis.plot(np.arange(1, Type.size+1), Avg_Disturbance * np.ones(Type.size), color='r', linestyle='--')
    #     Axis.text(0.5*Type.size+1, Avg_Disturbance, str(round((Avg_Disturbance-1)*100))+ '% Error', ha='right', va='bottom', color='r', fontsize=Settings.get('FontSize', 12))
    # # Legend
    # if Normalize:
    #     Axis.legend([Handle_Mean1, Handle_Mean1_Disturbance, Handle_Mean2_Disturbance], ['Without Disturbance', Legend[1], Legend[2]])
    # else:
    #     Axis.legend([Handle_Mean1, Handle_Mean1_Disturbance, Handle_Mean2, Handle_Mean2_Disturbance], Legend)
    # Fit Axis to Figure
    Figure.tight_layout()   
    # Return Handles
    Figure = Axis.figure
    return Figure, Axis

def Visualize_RPA2_Sorted(Threshold, Input_1, Input_2, Output_1, Output_1_Disturbance, Output_2, Output_2_Disturbance, Color, Color_Disturbance, Width, Title, Legend, Settings, Normalize, XLabels, Axis=None):
    # Create Figure and Axis
    if Axis is None:
        Figure, Axis = plt.subplots(figsize=Settings.get('FigureSize', (8, 4)))
    # Concatenate the Two Inputs/Outputs
    Input = pd.concat([Input_1, Input_2])
    Output_1.columns = Output_2.columns
    Output = pd.concat([Output_1, Output_2])
    Output_1_Disturbance.columns = Output_2_Disturbance.columns
    Output_Disturbance = pd.concat([Output_1_Disturbance, Output_2_Disturbance])
    Type = pd.DataFrame(np.concatenate([np.ones(Output_1.shape[0]), 2*np.ones(Output_2.shape[0])]))
    # Compute Means
    Mean = Output.mean(1)
    Mean_Disturbance = Output_Disturbance.mean(1)
    # Compute Standard Deviations
    Std = Output.std(1)
    Std_Disturbance = Output_Disturbance.std(1)
    # Filter
    Filtering_Flags = Mean <= Threshold
    Filtering_Flags = Filtering_Flags.values
    Input = Input.loc[Filtering_Flags]
    Output = Output.loc[Filtering_Flags]
    Output_Disturbance = Output_Disturbance.loc[Filtering_Flags]
    Mean = Mean.loc[Filtering_Flags]
    Mean_Disturbance = Mean_Disturbance.loc[Filtering_Flags]
    Std = Std.loc[Filtering_Flags]
    Std_Disturbance = Std_Disturbance.loc[Filtering_Flags]
    Type = Type.loc[Filtering_Flags]
    # Sort
    Sorting_Indeces = Mean.argsort()
    Input = Input.iloc[Sorting_Indeces]
    Output = Output.iloc[Sorting_Indeces]
    Output_Disturbance = Output_Disturbance.iloc[Sorting_Indeces]
    Mean = Mean.iloc[Sorting_Indeces]
    Mean_Disturbance = Mean_Disturbance.iloc[Sorting_Indeces]
    Std = Std.iloc[Sorting_Indeces]
    Std_Disturbance = Std_Disturbance.iloc[Sorting_Indeces]
    Type = Type.iloc[Sorting_Indeces]
    # Print Range
    print(f"Output Range: [{Mean.iloc[0]}, {Mean.iloc[-1]}],  Fold Change: {Mean.iloc[-1]/Mean.iloc[0]}.")
    # Normalize
    Normalize
    if Normalize:
        Output = Output.div(Mean, axis=0)
        Output_Disturbance = Output_Disturbance.div(Mean, axis=0)
        Output.div(Mean, axis=0)
        Std_Disturbance = Std_Disturbance / Mean
        Std = Std / Mean
        Mean_Disturbance = Mean_Disturbance / Mean
        Mean = Mean / Mean
    # Plot Means
    Sorted_Colors = np.where(Type == 1, Color[0], Color[1]) 
    Sorted_Colors_Disturbance = np.where(Type == 1, Color_Disturbance[0], Color_Disturbance[1]) 
    Handle_Mean1 = Axis.bar(np.arange(1, Input.size+1) - Width/2 - Width/10, Mean, Width, color=Sorted_Colors, edgecolor='black')
    Handle_Mean1_Disturbance = Axis.bar(np.arange(1, Input.size+1) + Width/2 + Width/10, Mean_Disturbance, Width, color=Sorted_Colors_Disturbance, edgecolor='black')
    # Error Bars
    Axis.errorbar(np.arange(1, Input.size+1) - Width/2 - Width/10, Mean, yerr=Std, fmt='None', capsize=4, ecolor='black', linestyle='None')
    Axis.errorbar(np.arange(1, Input.size+1) + Width/2 + Width/10, Mean_Disturbance, yerr=Std_Disturbance, fmt='None', capsize=4, ecolor='black', linestyle='None')
    # Plot Data Points
    Axis.plot(np.arange(1, Input.size+1) - Width/2 - Width/10, Output, marker='o', linestyle='', markeredgecolor='black', markerfacecolor=Color[0], markersize=Settings.get('MarkerSize', 2), markeredgewidth=0.2)        
    Axis.plot(np.arange(1, Input.size+1) + Width/2 + Width/10, Output_Disturbance, marker='o', linestyle='', markeredgecolor='black', markerfacecolor=Color_Disturbance[0], markersize=Settings.get('MarkerSize', 2), markeredgewidth=0.2)        
    # Axis Limits
    if 'XLimits' in Settings:
        Axis.set_xlim(Settings['XLimits'])
    if 'YLimits' in Settings:
        Axis.set_ylim(Settings['YLimits'])
    # Title
    Axis.set_title(Title, fontweight='bold')
    # Axis Labels
    Axis.set_xlabel('Arabinose (%)', fontsize = Settings.get('FontSize', 12), labelpad=70)
    if Normalize:
        Axis.set_ylabel('Normalized Output', fontsize = Settings.get('FontSize', 12))
    else:
        Axis.set_ylabel('Mean Output (MEF)', fontsize = Settings.get('FontSize', 12))
    # Axis Ticks
    Axis.set_xticks(np.arange(1, Input.size + 1))
    Axis.set_xticklabels([])
    heights = np.where(Type == 1, -Axis.get_ylim()[1]*0.2, -Axis.get_ylim()[1]*0.28)
    for i, height in enumerate(heights):
        Axis.text(i+1, height, str(Input.iloc[i]), ha='center', va='bottom', fontsize=Settings.get('FontSize', 12))
    heights_PlaceHolders = np.where(Type == 2, -Axis.get_ylim()[1]*0.17, -Axis.get_ylim()[1]*0.25)
    for i, height in enumerate(heights_PlaceHolders):
        Axis.text(i+1, height, '_', ha='center', va='bottom', fontsize=1.5*Settings.get('FontSize', 12))
    if not XLabels:
        Axis.set_xticklabels([])
    # Ticks Fontsize
    Axis.tick_params(axis='both', which='major', labelsize=Settings.get('FontSize', 12))
    # Grid
    Axis.set_axisbelow(True)
    Axis.grid(True)
    Axis.grid(color='gray', linestyle='--', linewidth=0.5)
    # Average Error
    if Normalize:
        Avg_Disturbance = Mean_Disturbance.mean()
        Axis.plot(np.arange(1, Input.size+1), Avg_Disturbance * np.ones(Input.size), color='r', linestyle='--')
        Axis.text(0.95*Input.size+1, Avg_Disturbance*1.095, str(round((Avg_Disturbance-1)*100))+ '% Error', ha='right', va='bottom', color='r', fontsize=Settings.get('FontSize', 12))
    # Legend
    # if Normalize:
    #     Axis.legend([Handle_Mean1, Handle_Mean1_Disturbance, Handle_Mean2_Disturbance], ['Without Disturbance', Legend[1], Legend[2]])
    # else:
    #     Axis.legend([Handle_Mean1, Handle_Mean1_Disturbance, Handle_Mean2, Handle_Mean2_Disturbance], Legend)
    # Fit Axis to Figure
    plt.tight_layout()  
    # Return Handles
    Figure = Axis.figure
    return Figure, Axis