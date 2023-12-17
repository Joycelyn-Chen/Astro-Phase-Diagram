# Astro-Phase-Diagram
- This is the notebook and the code where Joycelyn tried to generate the energy map for the magnetohydrodynamic synthesis dataset from [Alex Hill et al 2012]()

## `energy_map.ipynb`
- This is the experimental notebook where Joycelyn plays around the data and confirms which attribute to read, such as the density, temperature, velocity...etc. 
- Final code can be found in the bottom-most cell
![total energy map](./graphs/total_energy.png)
> Energy map for the entire slice

![masked enrygy map](./graphs/masked_energy.png)
> masked area of the energy map


## `phase_diagram.py`
- this is the complete code for Joyceyln to plot the phase diagram for each timestamp, ultimatly made into a movie
- It depicts the Thermal pressure as a function of density for each timestamp. (scatter chart)
- In order to get a better understanding of the region, we added the histogram graphs for density, temperature and the thermal pressure (all in log scale).

## `tracing_SN_events.ipynb`
- this is the organized notebook for tracing the origin of a single SN event. 
- The code allows you to input the desire region of interest and output all the SN that goes off within the region from the past 10 Myr