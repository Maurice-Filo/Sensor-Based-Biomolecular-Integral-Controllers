import FlowCytometryTools
import numpy as np
from FlowCytometryTools import FCMeasurement
import matplotlib.pyplot as plt

# Load the fcs files
Sample_OL = FCMeasurement(ID='TestSample', datafile='Nature_OL.fcs')
Sample_rAIF = FCMeasurement(ID='TestSample', datafile='Nature_rAIF.fcs')

# Colors
Color_OL = (0.5, 0.5, 0.5)
Color_rAIF = tuple(np.array([80, 137, 189]) / 255)
Color_sAIF = tuple(np.array([78, 146, 103]) / 255)
Opacity = 0.7

# Compute Statistics
Mean_OL = Sample_OL.data['GFP-A'].mean()
Mean_rAIF = Sample_rAIF.data['GFP-A'].mean()
CV_OL = Sample_OL.data['GFP-A'].std() / Mean_OL
CV_rAIF = Sample_rAIF.data['GFP-A'].std() / Mean_rAIF
Min_OL = Sample_OL.data['GFP-A'].min()
Min_rAIF = Sample_rAIF.data['GFP-A'].min()
Max_OL = Sample_OL.data['GFP-A'].max()
Max_rAIF = Sample_rAIF.data['GFP-A'].max()

# Bins
Bins = np.logspace(np.log10(10), np.log10(max(Max_OL, Max_rAIF)), num=1000)

# Figure for rAIF
Figure, Axis = plt.subplots(figsize=(8, 4))
# Histograms
plt.hist(Sample_OL.data['GFP-A'], bins=Bins, color=Color_OL, alpha=Opacity, label='OL')
plt.hist(Sample_rAIF.data['GFP-A'], bins=Bins, color=Color_rAIF, alpha=Opacity, label='rAIF')
# Axis Labels
plt.xlabel('Output (a.u.)', fontsize = 14)
plt.ylabel('Count', fontsize = 14)
# Axis Scale
plt.xscale('log')
# Axis Limits
plt.ylim([0, 600])
# Fontsize
plt.tick_params(axis='both', which='major', labelsize=14)
# Legend
plt.legend(loc='upper right', fontsize=14)
# Grid
Axis.set_axisbelow(True)
plt.grid(True)
plt.grid(color='gray', linestyle='--', linewidth=0.5)
# Means
plt.plot(Mean_OL * np.array([1,1]), [0, 600], linestyle='-', color=Color_OL)
plt.plot(Mean_rAIF * np.array([1,1]), [0, 600], linestyle='-', color=Color_rAIF)
# Annotations
plt.annotate('Mean = ' + str(round(Mean_OL)) + '\nCV = ' + str(round(CV_OL, 2)) + '\nCounts = ' + str(Sample_OL.data.shape[0]), 
             xy=(5e3, 200), xytext=(2e4, 300),
             arrowprops=dict(facecolor=Color_OL, edgecolor=Color_OL, shrink=0.05),
             horizontalalignment='left', fontsize=14)
plt.annotate('Mean = ' + str(round(Mean_rAIF)) + '\nCV = ' + str(round(CV_rAIF, 2)) + '\nCounts = ' + str(Sample_rAIF.data.shape[0]), 
             xy=(6e2, 200), xytext=(3e2, 300),
             arrowprops=dict(facecolor=Color_rAIF, edgecolor=Color_rAIF, shrink=0.05),
             horizontalalignment='right', fontsize=14)
# Fit Axis to Figure
plt.tight_layout() 
plt.show()
# Save Figure
Save_Flag = True
if Save_Flag:
    Figure.savefig("Published_rAIF.pdf")