import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Modules import Visualize_RPA1_Sorted
from Modules import Visualize_RPA2_Sorted

import sys
print(sys.executable)

# Specify the path to Excel file
file_path = 'Steady_state_error_triplicates_data.xlsx'

# Use pandas to read the Excel file
df = pd.read_excel(file_path)

# Extract Input Data (Arabinose)
Input = df.iloc[6:15, 1]

# Extract Open Loop Data
OL1_Weak = df.iloc[6:15, 2:5]
OL1_Weak_aTc = df.iloc[6:15, 5:8]
OL2_Weak = df.iloc[6:15, 8:11]
OL2_Weak_aTc = df.iloc[6:15, 11:14]
OL1_Strong = df.iloc[6:15, 14:17]
OL1_Strong_aTc = df.iloc[6:15, 17:20]
OL2_Strong = df.iloc[6:15, 20:23]
OL2_Strong_aTc = df.iloc[6:15, 23:26]

# Extract P Data
P_Weak = df.iloc[6:15, 26:29]
P_Weak_aTc = df.iloc[6:15, 29:32]
P_Strong = df.iloc[6:15, 38:41]
P_Strong_aTc = df.iloc[6:15, 41:44]

# Extract PI Data
PI_Weak = df.iloc[6:15, 32:35]
PI_Weak_aTc = df.iloc[6:15, 35:38]
PI_Strong = df.iloc[6:15, 44:47]
PI_Strong_aTc = df.iloc[6:15, 47:50]

# Plot Settings
Color_Undisturbed = tuple(np.array([[78, 146, 103, 255], [78, 146, 103, 100]]) / 255)
Color_Disturbed = tuple(np.array([[207, 92, 54, 255], [207, 92, 54, 100]]) / 255)
Width = 0.2
Opacity = 0.1
Plot_Settings = {
    'MarkerSize' : 3,
    'FontSize' : 16,
    'FigureSize' : (8,5),
}

# Open Loop
Input_OL = pd.DataFrame(np.concatenate([2*np.ones(2), np.ones(2)]))  # 1:Strong, 2:Weak
OL_Index = 0
OL1 = OL2_Weak.iloc[[OL_Index]]
OL2 = OL1_Weak.iloc[[OL_Index]]
OL3 = OL2_Strong.iloc[[OL_Index]]
OL4 = OL1_Strong.iloc[[OL_Index]]
OL2.columns = OL1.columns
OL3.columns = OL1.columns
OL4.columns = OL1.columns
OL1_aTc = OL2_Weak_aTc.iloc[[OL_Index]]
OL2_aTc = OL1_Weak_aTc.iloc[[OL_Index]]
OL3_aTc = OL2_Strong_aTc.iloc[[OL_Index]]
OL4_aTc = OL1_Strong_aTc.iloc[[OL_Index]]
OL2_aTc.columns = OL1_aTc.columns
OL3_aTc.columns = OL1_aTc.columns
OL4_aTc.columns = OL1_aTc.columns
OL = pd.concat([OL1, OL2, OL3, OL4], axis=0)
OL_aTc = pd.concat([OL1_aTc, OL2_aTc, OL3_aTc, OL4_aTc], axis=0)

Legend = ['Without Disturbance, Strong P$_{ARA}$', 'With Disturbance, Strong P$_{ARA}$', 'Without Disturbance, Weak P$_{ARA}$', 'With Disturbance, Weak P$_{ARA}$']

# OL Sorted
Width = 0.1
Threshold = 55
Plot_Settings['FigureSize'] = (8,4)
Figure_OL_Sorted, Axis_OL_Sorted = Visualize_RPA1_Sorted(Threshold, Input_OL, OL, OL_aTc, Color_Undisturbed, Color_Disturbed, Width, Title='', Legend=Legend, Settings=Plot_Settings, Normalize=False, XLabels=False, Axis=None)

# P Filtered and Sorted
Width = 0.3
Plot_Settings['FigureSize'] = (8,5.2)
Plot_Settings['YLimits'] = (0,55)
Figure_P_Sorted, Axis_P_Sorted = Visualize_RPA2_Sorted(Threshold, Input, Input, P_Strong, P_Strong_aTc, P_Weak, P_Weak_aTc, Color_Undisturbed, Color_Disturbed, Width, Title='', Legend=Legend, Settings=Plot_Settings, Normalize=False, XLabels=True, Axis=None)

# PI Filtered and Sorted
Figure_PI_Sorted, Axis_PI_Sorted = Visualize_RPA2_Sorted(Threshold, Input, Input, PI_Strong, PI_Strong_aTc, PI_Weak, PI_Weak_aTc, Color_Undisturbed, Color_Disturbed, Width, Title='', Legend=Legend, Settings=Plot_Settings, Normalize=False, XLabels=True, Axis=None)
plt.show()

# Save the Figures to PDF Files
Save_Flag = False
if Save_Flag:
    Figure_OL_Sorted.savefig("RPA_OL_Sorted.pdf")
    Figure_P_Sorted.savefig("RPA_P_Sorted.pdf")
    Figure_PI_Sorted.savefig("RPA_PI_Sorted.pdf")