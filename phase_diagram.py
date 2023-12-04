import matplotlib.pyplot as plt
import yt
import numpy as np
import os
import glob
from astropy import units as u
import cv2 as cv
os.environ.pop("QT_QPA_PLATFORM_PLUGIN_PATH")


def plot_histogram(y, category, bins = 1000):
    hist, bins, _ = plt.hist(y, bins=bins)
    logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
    plt.hist(x, bins=logbins)
    plt.title(f"Histogram of {category} (masked region)")
    plt.xlabel(f"{category}")
    plt.ylabel("# of bins")
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True) 

    plt.savefig(os.path.join(graph_root, "histogram", category, f"phase_{timestamp}_no_zero_init.png"))
    plt.clf()


def volume_fraction(dens, temp, thres_1 = pow(10, 4.2), thres_2 = pow(10, 5.5)):
    # volume = 2.9379989445851796e+55 cm**3 / 1.0000000000000004*yt.units.pc**3
    unit_vol = 1.0000000000000004 * yt.units.pc**3
    cold_mass = 0
    cold_vol = 0
    warm_mass = 0
    warm_vol = 0
    transition_mass = 0
    transition_vol = 0

    print("Calculating volume fraction...\n")
    for i, value in enumerate(temp):
        if value < thres_1:
            cold_mass += dens[i] * unit_vol
            cold_vol += unit_vol
        elif value < thres_2:
            transition_mass += dens[i] * unit_vol
            transition_vol += unit_vol
        else:
            warm_mass += dens[i] * unit_vol
            warm_vol += unit_vol
    return cold_mass, cold_vol, transition_mass, transition_vol, warm_mass, warm_vol



def plot_line(x, cool_y, transition_y, warm_y, mass_or_vol):
    plt.plot(x, cool_y, label = f"cool {mass_or_vol}") 
    plt.plot(x, transition_y, label = f"transition {mass_or_vol}") 
    plt.plot(x, warm_y, label = f"warm {mass_or_vol}") 
    
    plt.legend() 
    plt.title(f"Total {mass_or_vol}")
    plt.xlabel("time (Myr)")
    plt.ylabel(f"integrating {mass_or_vol}")
    # plt.xscale('log')
    plt.yscale('log')
    plt.grid(True) 

    plt.savefig(os.path.join(graph_root, "mass_fraction", f"{mass_or_vol}.png"))
    plt.clf()


# Initialization

# datadir = "/Users/joycelynchen/Desktop/UBC/Research/Data/synthetic_data/synthetic_data"
datadir = "/home/joy0921/Desktop/Dataset/200_210/raw_data"
mask_root = "/home/joy0921/Desktop/Dataset/200_210/case_masks"
filename_prefix = "sn34_smd132_bx5_pe300_hdf5_plt_cnt_0"
graph_root = "./graphs"


# Loading the data

begin_time = 202
end_time = 209

temp_arr = []
dens_arr = []
P_th_arr = []
mass_vol_values = []
mass_vol_values_list = [[], [], [], [], [], []]

for timestamp in range(begin_time, end_time + 1, 1):
    # Init array
    temp_arr.clear()
    dens_arr.clear()
    P_th_arr.clear()


    print(f"Processing timestamp {timestamp}")
    ds = yt.load(os.path.join(datadir, f"{filename_prefix}{timestamp}"))

    center =  [0, 0, 0]*yt.units.pc
    arb_center = ds.arr(center,'code_length')
    left_edge = arb_center - ds.quan(400,'pc')
    right_edge = arb_center + ds.quan(400,'pc')
    obj = ds.arbitrary_grid(left_edge, right_edge, dims=[800, 800, 800]*yt.units.pc)



    # Initializing value storing entity
    single_phase = yt.YTArray(np.zeros((800, 800, 15))) * (yt.units.kelvin) / yt.units.pc**3
    dens = yt.YTArray(np.zeros((800, 800, 15))) * (1) / yt.units.pc**3
    masked_single_phase = yt.YTArray(np.zeros((800, 800, 15))) * (yt.units.kelvin) / yt.units.pc**3
    masked_dens = yt.YTArray(np.zeros((800, 800, 15))) * (1) / yt.units.pc**3
    masked_temp = yt.YTArray(np.zeros((800, 800, 15))) * (yt.units.kelvin)



    # Thermal pressure = n * temp
    mu = 1.4
    m_H = yt.physical_constants.mass_hydrogen


    middle_plane = 300
    start = middle_plane - 7
    end = middle_plane + 8


    for k in range(start, end, 1):
        temp = obj["flash", "temp"][:, :, k]
        n = obj["flash", "dens"][:,:,k] / (mu * m_H) 
        dens[:, :, k-start] = (n).to(dens.units)
        single_phase[:, :, k-start] = (n * temp).to(single_phase.units)
        print(f"Whole slice {k} done. [ {round(((k-start)/15)*100, 1)}% ]")

    #single_sum = single_phase.sum(axis=2)
    #dens_sum = dens.sum(axis=2)

    for z in range(start, end, 1):
        mask_filename = f"{filename_prefix}{timestamp}_z{z}.png"
        mask = cv.imread(os.path.join(mask_root, str(timestamp), mask_filename))
        coordinates = np.argwhere(mask == 255)
        for coord in coordinates:
            x, y = coord[0], coord[1]
            temp = obj["flash", "temp"][x,y,z]
            n = obj["flash", "dens"][x,y,z] / (mu * m_H)
            masked_dens[x, y, z-start] = (n).to(masked_dens.units)
            masked_temp[x, y, z-start] = temp.to(masked_temp.units)
            masked_single_phase[x, y, z-start] = (n * temp).to(masked_single_phase.units)
            temp_arr.append(temp)
            dens_arr.append(n)
            P_th_arr.append(n*temp)

            
        print(f"Masked region slice {z} done. [ {round(((z-start)/15)*100, 1)}% ]")

    # masked_single_sum = masked_single_phase.sum(axis=2)
    # masked_dens_sum = masked_dens.sum(axis=2)


    # thermal pressure as a function of density in scatter plot
    dens_flattened = np.array(dens.to_ndarray().flatten())
    P_th_flattened = np.array(single_phase.to_ndarray().flatten())
    masked_dens_flattened = np.array(masked_dens.to_ndarray().flatten())
    masked_P_th_flattened = np.array(masked_single_phase.to_ndarray().flatten())
    masked_temp_flattened = np.array(masked_temp.to_ndarray().flatten())

    cold_mass, cold_vol, transition_mass, transition_vol, warm_mass, warm_vol = volume_fraction(masked_dens_flattened, masked_temp_flattened, thres_1 = 10**4.2, thres_2 = 10**5.5)
    mass_vol_values.append((cold_mass, cold_vol, transition_mass, transition_vol, warm_mass, warm_vol))

    mass_vol_values_list[0].append(cold_mass)
    mass_vol_values_list[1].append(cold_vol)
    mass_vol_values_list[2].append(transition_mass)
    mass_vol_values_list[3].append(transition_vol)
    mass_vol_values_list[4].append(warm_mass)
    mass_vol_values_list[5].append(warm_vol)
    
    print(f"Mass and volume: ")
    print(f"Cold mass: {mass_vol_values_list[0]}")
    print(f"Cold volume: {mass_vol_values_list[1]}")
    print(f"transition mass: {mass_vol_values_list[2]}")
    print(f"transition volume: {mass_vol_values_list[3]}")
    print(f"Warm mass: {mass_vol_values_list[4]}")
    print(f"Warm volume: {mass_vol_values_list[5]}")

    # Scatter plot for dens - thermal pressure
    # color: gray(#7f7f7f), 
    # plt.scatter(x, y, c='#1f77b4', marker='o', label='whole slice' if timestamp == begin_time else "", s=0.01)
    # plt.scatter(x, y1, c='r', marker='o', label='masked area' if timestamp == begin_time else "", s=0.01)
    
    # # Plotting scatters
    # plt.scatter(dens_flattened, P_th_flattened, c='#1f77b4', marker='o', label='whole slice', s=0.01)
    # plt.scatter(dens_flattened, masked_P_th_flattened, c='r', marker='o', label='masked area', s=0.01)
    # plt.legend(loc='upper left')

    # plt.title("Single Phase Diagram")
    # plt.xlabel("Density (n)")
    # plt.ylabel("Thermal pressure (nT)")
    # plt.grid(True) 
    # plt.yscale('log')
    # plt.xscale('log')  
    # # plt.show() 
    # plt.savefig(os.path.join(graph_root, f"phase_{timestamp}_1.png"))
    # plt.clf()


    # # histogram of thermal pressure
    # plot_histogram(P_th_arr, "thermal_pressure", bins = 1000)

    # # histogram of density
    # plot_histogram(dens_arr, "density", bins = 1000)

    # # histogram of temperature
    # plot_histogram(temp_arr, "temperature", bins = 1000)


# for item in mass_vol_values:
#     mass_vol_values_list[0].append(item[0])
#     mass_vol_values_list[1].append(item[1])
#     mass_vol_values_list[2].append(item[2])
#     mass_vol_values_list[3].append(item[3])
#     mass_vol_values_list[4].append(item[4])
#     mass_vol_values_list[5].append(item[5])
# mass_vol_values_list = map(list, zip(*mass_vol_values))




x_range = range(begin_time, end_time + 1)

# plot lines 
plot_line(x_range, mass_vol_values_list[0], mass_vol_values_list[2], mass_vol_values_list[4], 'mass')
plot_line(x_range, mass_vol_values_list[1], mass_vol_values_list[3], mass_vol_values_list[5], 'volume')
