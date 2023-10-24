import matplotlib.pyplot as plt
import yt
import numpy as np
import os
import glob
from astropy import units as u
import cv2 as cv
os.environ.pop("QT_QPA_PLATFORM_PLUGIN_PATH")

# Initialization

# datadir = "/Users/joycelynchen/Desktop/UBC/Research/Data/synthetic_data/synthetic_data"
datadir = "/home/joy0921/Desktop/2023S/Dataset/200_360/finer_time_200_360_original"
mask_root = "/home/joy0921/Desktop/2023S/Dataset/200_210/case_masks"
filename_prefix = "sn34_smd132_bx5_pe300_hdf5_plt_cnt_0"
graph_root = "./graphs"


# Loading the data

begin_time = 201
end_time = 203


for timestamp in range(begin_time, end_time, 1):
    print(f"Processing timestamp {timestamp}")
    ds = yt.load(os.path.join(datadir, f"{filename_prefix}{timestamp}"))

    center =  [0, 0, 0]*yt.units.pc
    arb_center = ds.arr(center,'code_length')
    left_edge = arb_center - ds.quan(400,'pc')
    right_edge = arb_center + ds.quan(400,'pc')
    obj = ds.arbitrary_grid(left_edge, right_edge, dims=[800, 800, 800]*yt.units.pc)



    # Initializing value storing entity
    single_phase = yt.YTArray(np.zeros((800, 800, 15))) * (yt.units.kelvin) / yt.units.cm**3
    dens = yt.YTArray(np.zeros((800, 800, 15))) * (1) / yt.units.cm**3
    masked_single_phase = yt.YTArray(np.zeros((800, 800, 15))) * (yt.units.kelvin) / yt.units.cm**3
    masked_dens = yt.YTArray(np.zeros((800, 800, 15))) * (1) / yt.units.cm**3



    # Thermal pressure = n * temp
    mu = 1.4
    m_H = yt.physical_constants.mass_hydrogen


    middle_plane = 300
    start = middle_plane - 7
    end = middle_plane + 8

    for k in range(start, end, 1):
        for i in range(800):
            for j in range(800):
                temp = obj["flash", "temp"][i,j,k]
                n = obj["flash", "dens"][i,j,k] / (mu * m_H)
                dens[i][j][k-start] = (n).to(dens.units)
                single_phase[i][j][k-start] = (n * temp).to(single_phase.units)
        print(f"Whole slice {k} done. [ {round(((k-start)/15)*100, 1)}% ]")
        
    single_sum = single_phase.sum(axis=2)
    dens_sum = dens.sum(axis=2)

    for z in range(start, end, 1):
        mask_filename = f"{filename_prefix}{timestamp}_z{z}.png"
        mask = cv.imread(os.path.join(mask_root, str(timestamp), mask_filename))
        coordinates = np.argwhere(mask == 255)
        for coord in coordinates:
            x, y = coord[0], coord[1]
            temp = obj["flash", "temp"][x,y,z]
            n = obj["flash", "dens"][x,y,z] / (mu * m_H)
            masked_dens[i][j][z-start] = (n).to(masked_dens.units)
            masked_single_phase[x][y][z-start] = (n * temp).to(masked_single_phase.units)

            
        print(f"Masked region slice {z} done. [ {round(((z-start)/15)*100, 1)}% ]")

    masked_single_sum = masked_single_phase.sum(axis=2)
    masked_dens_sum = masked_dens.sum(axis=2)


    # thermal pressure as a function of density in scatter plot
    x = np.array(dens_sum.to_ndarray().flatten())
    y = np.array(single_sum.to_ndarray().flatten())
    x1 = np.array(masked_dens_sum.to_ndarray().flatten())
    y1 = np.array(masked_single_sum.to_ndarray().flatten())


    # color: gray(#7f7f7f), 
    plt.scatter(x, y, c='#1f77b4', marker='o', label='whole slice' if timestamp == begin_time else "", s=0.01)
    plt.scatter(x, y1, c='r', marker='o', label='masked area' if timestamp == begin_time else "", s=0.01)
    plt.legend(loc='upper left')

    plt.title("Single Phase Diagram")
    plt.xlabel("Density (n)")
    plt.ylabel("Thermal pressure (nT)")
    plt.grid(True) 
    plt.yscale('log')
    plt.xscale('log')  
    # plt.show() 
    plt.savefig(os.path.join(graph_root, f"phase_{timestamp}_tmp.png"))
