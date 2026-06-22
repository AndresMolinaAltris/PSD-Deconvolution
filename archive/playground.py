import matplotlib.pyplot as plt
import numpy as np


#####################################################################
## SMALL PARTICLE POPULATION STATISTICS #############################
##################################################################### 
y1_50 = [2.64, 2.122, 2.645, 2.608, 2.331, 2.221, 2.644, 2.645, np.nan, np.nan, np.nan,
      2.765, 1.824, 2.469, 2.249, 2.283, 2.327, 2.838, 2.618]
y1_10 = [1.901, 1.631, 1.907, 1.894, 1.757, 1.68, 1.915, 1.907, np.nan, np.nan, np.nan,
        1.954, 1.212, 1.752, 1.602, 1.63, 1.555, 1.553, 1.604]

y1_90 = [3.666, 2.76, 3.669, 3.59, 3.095, 2.937, 3.65, 3.669, np.nan, np.nan, np.nan,
        3.913, 2.743, 3.478, 3.159, 3.199, 3.484, 5.186, 4.274]


#####################################################################
## BIG PARTICLE POPULATION STATISTICS #############################
##################################################################### 
y2_50 = [27.106, 23.185, 26.102, 25.049, 25.397, 23.497, 27.543, 26.17,
        22.447, 20.387, 23.102, 27.046, 16.878, 19.609, 20.495, 19.481,
        19.308, 22.063, 22.82]

y2_10 = [14.091, 12.149, 13.402, 12.932, 12.981, 12.183, 14.02, 13.439, 11.526, 10.357,
        11.498, 13.647, 8.583, 9.64, 10.44, 9.547, 9.378, 10.818, 11.255]

y2_90 = [52.143, 44.242, 50.836, 48.522, 49.689, 45.319, 54.108, 50.961, 43.714, 40.13,
        46.415, 53.601, 33.19, 39.89, 40.232, 39.751, 39.753, 44.998, 46.268]



# X values from 1 to len(y)
x = list(range(1, len(y1) + 1))

###########################################################
###########################################################
#### PLOT 1 ###############################################
plt.scatter(x, y1_90, facecolors='none', edgecolors='blue', marker='D', label=r'$D_{90}$')
plt.scatter(x, y1_50, color='blue', label=r'$D_{50}$')
plt.scatter(x, y1_10, facecolors='none', edgecolors='blue', label=r'$D_{10}$')
plt.title("Small Particle Population")
plt.xlabel("Sample Number")
plt.xticks(range(1, len(y1) + 1))  # Set ticks to integers
plt.ylabel(r"Particle Diameter / $\mu$m")
plt.legend()
plt.grid(True, color='gray', linestyle='--', alpha=0.5)
plt.show()
#### PLOT 2 ###############################################
plt.scatter(x, y2_90, facecolors='none', edgecolors='green', marker='D', label=r'$D_{90}$')
plt.scatter(x, y2_50, color='green', label=r'$D_{50}$')
plt.scatter(x, y2_10, facecolors='none', edgecolors='green', label=r'$D_{10}$')
plt.title("Large Particle Population")
plt.xlabel("Sample Number")
plt.ylabel(r"Particle Diameter / $\mu$m")
plt.xticks(range(1, len(y2_90) + 1))  # Set ticks to integers
plt.legend()
plt.grid(True, color='gray', linestyle='--', alpha=0.5)
plt.show()
###########################################################
###########################################################
###########################################################


###########################################################
# POPULATION PERCENTAGE ###################################
###########################################################
population_1 = [0.8, 0.4, 1.0, 1.0, 0.4, 0.4, 0.8, 1.0, 0, 0, 0,
                1.5, 1.3, 1.6, 1.3, 1.6, 2.2, 7.1, 4.4]

population_2 = [99.2, 99.6, 99.0, 99.0, 99.6, 99.6, 99.2, 99.0, 100, 100, 100,
                98.5, 98.6, 98.4, 98.7, 98.4, 97.8, 92.9, 95.6]


# X values from 1 to len(y)
x_pop = list(range(1, len(population_1) + 1))
#### PLOT 1 ###############################################
plt.scatter(x_pop, population_1, color='blue')
plt.title("Small Particle Population Percentage")
plt.xlabel("Sample Number")
plt.xticks(range(1, len(population_1) + 1))  # Set ticks to integers
plt.ylabel('Percentage / %')
plt.grid(True, color='gray', linestyle='--', alpha=0.5)
plt.show()
#### PLOT 2 ###############################################
plt.scatter(x_pop, population_2, color='green')
plt.title("Large Particle Population Percentage")
plt.xlabel("Sample Number")
plt.ylabel('Percentage / %')
plt.xticks(range(1, len(population_1) + 1))  # Set ticks to integers
plt.grid(True, color='gray', linestyle='--', alpha=0.5)
plt.show()
###########################################################
###########################################################
###########################################################

bet_data = [0.57, 0.46, 0.48, 0.59, 0.69, 0.5, 0.56, 0.577, 0.62,
            0.68, 0.72, 0.6, 0.89, 0.83, 1.33, 1.09]

population_1 = [0.8, 0.4, 1.0, 1.0, 0.4, 0.4, 0.8, 1.0, 
                1.5, 1.3, 1.6, 1.3, 1.6, 2.2, 7.1, 4.4]

plt.scatter(population_1, bet_data, color='blue')
plt.title("BET Correlation")
plt.xlabel("Population Percentage / %")
#plt.xticks(range(1, len(population_1) + 1))  # Set ticks to integers
plt.ylabel('BET SUrface Area (m2/g)')
plt.grid(True, color='gray', linestyle='--', alpha=0.5)
plt.show()

