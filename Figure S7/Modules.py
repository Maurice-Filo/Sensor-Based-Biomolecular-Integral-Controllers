import pandas as pd
import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
import pdb

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