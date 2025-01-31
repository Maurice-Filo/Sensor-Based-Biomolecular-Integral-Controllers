import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import trim_and_af_correct as trimaf

# White cell specification
fc_data_path_w = ['20230714_FlowCal_output.xlsx', '20230816_FlowCal_output.xlsx', '20230818_FlowCal_output.xlsx']
sample_ids_w = ['S001']
channel_w = 'FL4-H' 

# Sample specification
fc_data_path = ['20230714_FlowCal_output.xlsx', '20230816_FlowCal_output.xlsx', '20230818_FlowCal_output.xlsx']
channel = 'FL4-H' 
sample_ids = [[['S117'], 3], [['S134'], 3], [['S131'], 1]]
labels = ['Open Loop', 'P-Control', 'sAIF Control']

# Trimming algorithm parameters
## White cells
thresh_w = 0.005
bw_w = 0.05
## Samples
thresh = 0.005
bw = 0.05

# Trim & correct auto-fluoresence 
data_w = [None] * len(sample_ids)
trimmed_w = [None] * len(sample_ids)
data = [None] * len(sample_ids)
trimmed = [None] * len(sample_ids)
af_corrected = [None] * len(sample_ids)
means = [None] * len(sample_ids)
cvs = [None] * len(sample_ids)
means_trimmed = [None] * len(sample_ids)
cvs_trimmed = [None] * len(sample_ids)
means_w = [None] * len(sample_ids)
for i in range(len(sample_ids)):
    data_w[i] = trimaf.import_hists(fc_data_path_w[sample_ids[i][1]-1], sample_ids_w, channel_w)
    data[i] = trimaf.import_hists(fc_data_path[sample_ids[i][1]-1], sample_ids[i][0], channel)
    trimmed_w[i] = trimaf.trim(data_w[i], thresh=thresh_w, bw=bw_w)
    trimmed[i] = trimaf.trim(data[i], thresh=thresh, bw=bw)
    af_corrected[i] = trimaf.af_correct(trimmed[i], trimmed_w[i])
    means[i] = np.mean(data[i].events)
    cvs[i] = np.std(data[i].events) / np.mean(data[i].events)
    means_trimmed[i] = np.mean(trimmed[i].events)
    cvs_trimmed[i] = np.std(trimmed[i].events) / np.mean(trimmed[i].events)
    means_w[i] = np.mean(data_w[i].events)

# Plotting parameters
Opacity = [0.6, 0.4, 0.3]
Color_OL = np.array([244, 1, 154]) / 255
Color_sAIF = tuple(np.array([80, 137, 189]) / 255)
Color_P = tuple(np.array([255, 127, 80]) / 255)
Colors = [Color_OL, Color_P, Color_sAIF]

# Figure for data
Figure, Axis = plt.subplots(figsize=(10, 4))
# Auto-fluorescence
plt.hist(data_w[2].events[0], data_w[2].bins[0], edgecolor = (0, 0, 0), label='Autofluorescence', histtype='step')
plt.plot(means_w[2] * np.array([1,1]), [0, 550], linestyle='-', color=(0, 0, 0))
# Fluorescence
for i in range(len(sample_ids)):
    plt.hist(data[i].events[0], data[i].bins[0], color=Colors[i], alpha=Opacity[i], label=labels[i])
for i in range(len(sample_ids)):
    plt.hist(data[i].events[0], data[i].bins[0], edgecolor=Colors[i], histtype='step')
    plt.plot(means[i] * np.array([1,1]), [0, 550], linestyle='-', color=Colors[i])
# Axis Labels
plt.xlabel('Output (MEF)', fontsize = 14)
plt.ylabel('Count', fontsize = 14)
# Axis Scale
plt.xscale('log')
# Axis Limits
plt.ylim([0, 550])
plt.xlim([3,65])
# Axis Ticks
plt.xticks([5, 10, 25, 55], ['5', '10', '25', '55'])
# Fontsize
plt.tick_params(axis='both', which='major', labelsize=14)
# Legend
plt.legend(loc='upper right', fontsize=14)
# Grid
Axis.set_axisbelow(True)
plt.grid(True)
plt.grid(color='gray', linestyle='--', linewidth=0.5)
# Fit Axis to Figure
plt.tight_layout() 
plt.show()

# Save the Figures to PDF Files
Save_Flag = False
if Save_Flag:
    Figure.savefig("Distributions.pdf")

