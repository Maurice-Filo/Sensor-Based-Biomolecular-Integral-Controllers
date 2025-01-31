import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Modules import Visualize_Dynamics

# Specify the path to Excel file
file_path = 'Dynamics_triplicates_data.xlsx'

# Use pandas to read the Excel file
df = pd.read_excel(file_path)

# Extract Time Horizon
time = df.iloc[1:, 0]

# Extract Open Loop Data
OL = df.iloc[1:, 1:4]
OL_aTc = df.iloc[1:, 4:7]

# Extract P Data
P_High = df.iloc[1:, 7:10]
P_High_aTc = df.iloc[1:, 10:13]
P_Low = df.iloc[1:, 13:16]
P_Low_aTc = df.iloc[1:, 16:19]

# Extract PI Data
PI = df.iloc[1:, 19:22]
PI_aTc = df.iloc[1:, 22:25]

# Plot Settings
Color_Undisturbed = tuple(np.array([79, 146, 103]) / 255)
Color_Disturbed = tuple(np.array([207, 92, 54]) / 255)
Opacity = 0.1
Plot_Settings = {
    'MarkerSize' : 3,
    'XLimits' : [0, 7.1],
    'YLimits' : [0, 35], 
    'FigureSize' : (6.4, 4), 
    'FontSize' : 14,
    'Color_Disturbance' : tuple(np.array([241, 187, 97]) / 255)
}

# Open Loop
Figure_OL, Axis_OL = Visualize_Dynamics(time, OL, OL_aTc, Color_Undisturbed, Color_Disturbed, Opacity, '', Plot_Settings)

# P
Figure_P, Axis_P = Visualize_Dynamics(time, P_High, P_High_aTc, Color_Undisturbed, Color_Disturbed, Opacity, '', Plot_Settings)
Figure_P, Axis_P = Visualize_Dynamics(time, P_Low, P_Low_aTc, Color_Undisturbed, Color_Disturbed, Opacity, '', Plot_Settings, Axis_P)

# PI
Figure_PI, Axis_PI = Visualize_Dynamics(time, PI, PI_aTc, Color_Undisturbed, Color_Disturbed, Opacity, '', Plot_Settings)
plt.show()

# Save the Figures to PDF Files
Save_Flag = True
if Save_Flag:
    Figure_OL.savefig("Dynamics_OL.pdf")
    Figure_P.savefig("Dynamics_P.pdf")
    Figure_PI.savefig("Dynamics_PI.pdf")

print(Axis_OL.xaxis.label.get_fontname())