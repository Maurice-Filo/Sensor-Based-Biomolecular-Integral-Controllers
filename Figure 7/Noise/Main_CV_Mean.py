import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Modules import Visualize_Dynamics
from Modules import Visualize_CV_Mean

# Specify the path to Excel file
file_path = 'Mean_vs_CV_triplicates_data.xlsx'

# Use pandas to read the Excel file
df = pd.read_excel(file_path)

# Extract Input Data (Arabinose)
Input = df.iloc[6:15, 1]

# Extract Open Loop Data
OL2_Strong_Mean = df.iloc[6:15, 2:5]
OL2_Strong_CV = df.iloc[6:15, 5:8]
OL2_Weak_Mean = df.iloc[6:15, 8:11]
OL2_Weak_CV = df.iloc[6:15, 11:14]
OL1_Strong_Mean = df.iloc[6:15, 14:17]
OL1_Strong_CV = df.iloc[6:15, 17:20]
OL1_Weak_Mean = df.iloc[6:15, 20:23]
OL1_Weak_CV = df.iloc[6:15, 23:26]

# Extract P Data
P_Weak_Mean = df.iloc[6:15, 26:29]
P_Weak_CV = df.iloc[6:15, 29:32]
P_Strong_Mean = df.iloc[6:15, 38:41]
P_Strong_CV = df.iloc[6:15, 41:44]

# Extract PI Data
PI_Weak_Mean = df.iloc[6:15, 32:35]
PI_Weak_CV = df.iloc[6:15, 35:38]
PI_Strong_Mean = df.iloc[6:15, 44:47]
PI_Strong_CV = df.iloc[6:15, 47:50]

# Plot Settings
# Color_OL = (0.5, 0.5, 0.5)
# Color_P = tuple(np.array([207, 92, 54]) / 255)
# Color_PI = tuple(np.array([78, 146, 103]) / 255)
Color_OL = np.array([244, 1, 154]) / 255
Color_PI = tuple(np.array([80, 137, 189]) / 255)
Color_P = tuple(np.array([255, 127, 80]) / 255)
# Color_P = tuple(np.array([78, 146, 103]) / 255)
Opacity = 1
Opacity_High = 1
Opacity_MediumHigh= 0.5
Opacity_MediumLow= 0.3
Opacity_Low = 0.1
Plot_Settings = {
    'MarkerSize' : 2,
    'XLimits' : [0, 100],
    'YLimits' : [0, 2], 
    'FontSize' : 14,
    'FigureSize' : (10, 4),
}

## With Statistics
# Open Loop
Figure_All_Stats, Axis_All_Stats, Points_All_Stats = Visualize_CV_Mean(OL2_Strong_Mean, OL2_Strong_CV, Color_OL, Opacity_MediumHigh, 'Noise', Plot_Settings, True, Axis=None)
Visualize_CV_Mean(OL2_Weak_Mean, OL2_Weak_CV, Color_OL, Opacity_High, 'Noise', Plot_Settings,  True, Axis=Axis_All_Stats)
Visualize_CV_Mean(OL1_Weak_Mean, OL1_Weak_CV, Color_OL, Opacity_MediumLow, 'Noise', Plot_Settings, True, Axis=Axis_All_Stats)
Visualize_CV_Mean(OL1_Strong_Mean, OL1_Strong_CV, Color_OL, Opacity_Low, 'Noise', Plot_Settings, True,  Axis=Axis_All_Stats)

# P
Visualize_CV_Mean(P_Weak_Mean, P_Weak_CV, Color_P, Opacity_High, 'Noise', Plot_Settings, True,  Axis=Axis_All_Stats)
Visualize_CV_Mean(P_Strong_Mean, P_Strong_CV, Color_P, Opacity_MediumLow, 'Noise', Plot_Settings, True,  Axis=Axis_All_Stats)

# PI
Visualize_CV_Mean(PI_Weak_Mean, PI_Weak_CV, Color_PI, Opacity_High, 'Noise', Plot_Settings, True,  Axis=Axis_All_Stats)
Visualize_CV_Mean(PI_Strong_Mean, PI_Strong_CV, Color_PI, Opacity_MediumLow, 'Noise', Plot_Settings, True,  Axis=Axis_All_Stats)

## Without Statistics
Plot_Settings['MarkerSize'] = 7

# Open Loop
Opacity = 1
Figure_All, Axis_All, Points_OL2_Strong = Visualize_CV_Mean(OL2_Strong_Mean, OL2_Strong_CV, Color_OL, Opacity, '', Plot_Settings, False, Axis=None)
Figure_All, Axis_All, Points_OL2_Weak = Visualize_CV_Mean(OL2_Weak_Mean, OL2_Weak_CV, Color_OL, Opacity, '', Plot_Settings,  False, Axis=Axis_All)
Figure_All, Axis_All, Points_OL1_Weak = Visualize_CV_Mean(OL1_Weak_Mean, OL1_Weak_CV, Color_OL, Opacity, '', Plot_Settings, False, Axis=Axis_All)
Figure_All, Axis_All, Points_OL1_Strong  = Visualize_CV_Mean(OL1_Strong_Mean, OL1_Strong_CV, Color_OL, Opacity, '', Plot_Settings, False,  Axis=Axis_All)

# P
Figure_All, Axis_All, Points_P_Weak = Visualize_CV_Mean(P_Weak_Mean, P_Weak_CV, Color_P, Opacity, '', Plot_Settings, False,  Axis=Axis_All)
Figure_All, Axis_All, Points_P_Strong = Visualize_CV_Mean(P_Strong_Mean, P_Strong_CV, Color_P, Opacity, '', Plot_Settings, False,  Axis=Axis_All)

# PI
Figure_All, Axis_All, Points_PI_Weak = Visualize_CV_Mean(PI_Weak_Mean, PI_Weak_CV, Color_PI, Opacity, '', Plot_Settings, False,  Axis=Axis_All)
Figure_All, Axis_All, Points_PI_Strong = Visualize_CV_Mean(PI_Strong_Mean, PI_Strong_CV, Color_PI, Opacity, '', Plot_Settings, False,  Axis=Axis_All)
Axis_All.legend([Points_OL2_Strong[0], Points_P_Weak[0], Points_PI_Weak[0]], 
                ['Open Loop', 'P-Control', 'sAIF Control'], prop={'size': Plot_Settings.get('FontSize', 12)})
Axis_All.set_xlim([0,55])

# Function to be called when the mouse is moved
def on_move(event):
    if event.inaxes == Axis_All:
        Axis_All.set_title(f'x={event.xdata:.2f}, y={event.ydata:.2f}')
        Figure_All.canvas.draw()

# Connect the event to the function
Figure_All.canvas.mpl_connect('motion_notify_event', on_move)

plt.show()

# Save the Figures to PDF Files
Save_Flag = False
if Save_Flag:
    Figure_All.savefig("CV_Mean.pdf")