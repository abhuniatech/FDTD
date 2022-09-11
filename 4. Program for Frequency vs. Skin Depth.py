#!/usr/bin/env python
# coding: utf-8

# In[5]:


from math import pi, sqrt, atan
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

number = 7
number_of_frequencies = 8
permeability = 4*pi*10e-7

epsr = np.zeros((number_of_frequencies,number))
sigma = np.ones((number_of_frequencies,number))

angle = np.zeros(number_of_frequencies)
reflectance = np.zeros(number_of_frequencies)
skin_depth = np.zeros(number_of_frequencies)
n = np.zeros(number_of_frequencies)

epsr = [[1,72.9,12.7,66.0,55.8,31.6,80.1],   # 100 MHz
        [1,55.7,12.0,60.2,49.2,26.6,65.1],   # 200 MHz
        [1,49.8,11.7,58.2,46.8,24.8,60.0],   # 300 MHz
        [1,46.8,11.6,57.1,45.5,23.8,57.4],   # 400 MHz
        [1,44.9,11.5,56.4,44.6,23.2,55.8],   # 500 MHz
        [1,43.6,11.5,56.0,44.0,22.8,54.7],   # 600 MHz
        [1,42.7,11.4,55.6,43.5,22.5,53.9],   # 700 MHz
        [1,42.0,11.4,55.3,43.0,22.2,53.3]]   # 800 MHz

sigma = [[0.001,0.491,0.0684,0.708,0.475,0.306,0.559],
         [0.001,0.582,0.0726,0.743,0.518,0.335,0.639],
         [0.001,0.641,0.0765,0.771,0.553,0.356,0.692],
         [0.001,0.688,0.0807,0.796,0.586,0.374,0.737],
         [0.001,0.728,0.0854,0.822,0.621,0.391,0.779],
         [0.001,0.765,0.0905,0.850,0.658,0.407,0.819],
         [0.001,0.800,0.0962,0.879,0.697,0.423,0.860],
         [0.001,0.834,0.102,0.910,0.738,0.440,0.900]]

freq = [100,200,300,400,500,600,700,800]

print("Skin Depth vs. Frequency for different tissues")
print(" ")
print("Tissue 0 : Air,Tissue 1 : Skin,Tissue 2 : Fat,Tissue 3 :\
       Muscle,Tissue 4 : Cartilage,Tissue 5 : Lung,Tissue 0 : Brain ")
print(" ")

for j in range(0,7):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_title('Tissue %d' %j )
    ax.set(xlabel="Frequency",ylabel="Skin Depth")
    for k in range(0,8):
        omega = 2*pi*freq[k]
        skin_depth[k] = sqrt(2/(sigma[k][j]*permeability*omega))
    ax.plot(freq,skin_depth)
    plt.show()


# In[ ]:





# In[ ]:
