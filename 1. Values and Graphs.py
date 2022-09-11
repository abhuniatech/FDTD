#!/usr/bin/env python
# coding: utf-8

# In[105]:


from math import pi, sqrt, atan
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

number = 7
number_of_frequencies = 8
permeability = 4*pi*10e-7

epsr = np.zeros((number_of_frequencies,number))
sigma = np.ones((number_of_frequencies,number))

angle = np.zeros(number)
reflectance = np.zeros(number)
skin_depth = np.zeros(number)
n = np.zeros(number)

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

for j in range(number_of_frequencies):
    x=[]
    y=[]
    z=[]
    omega = 2*pi*freq[j]
    print("Frequency = ",freq[j],"MHz")
    for k in range(number):
        n[k] = sqrt(epsr[j][k])

        angle[k] = atan(n[k])*(180/pi)
        reflectance[k] = ((1-n[k])/(1+n[k]))**2
        skin_depth[k] = sqrt(2/(sigma[j][k]*permeability*omega))

        x.append(angle[k])
        y.append(reflectance[k])
        z.append(skin_depth[k])

    for k in range(number):
        df = pd.DataFrame({"Tissue":["Air","Skin","Fat","Muscle",
                           "Cartilage","Lung","Brain(GM)"],
                           "Brewster Angle":x,
                           "Reflectance":y,
                           "Skin Depth":z})


    display(df)
    print("    ")

for j in range(number):
    print("Frequency = ",freq[j],"MHz")
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set(xlabel="Tissues",ylabel="Skin Depth")
    skd = []
    omega = 2*pi*freq[j]
    for k in range(number):
        skd.append(sqrt(2/(sigma[j][k]*permeability*omega)))
        langs = ["Air","Skin","Fat","Muscle","Cartilage","Lung","Brain(GM)"]
    ax.bar(langs,skd)
    plt.show()



# In[90]:


print ("Tabulated Values ::")
print (" ")
print ("Frequency = 100 MHz")

df = pd.DataFrame({"Tissue":["Air","Skin","Fat",
                             "Muscle","Cartilage","Lung","Brain(GM)"],
                             "Epsilon(Permittivity)"
                             : [1,72.9,12.7,66.0,55.8,31.6,80.1],
                             "Sigma(Conductivity)"
                             : [0.001,0.491,0.0684,0.708,0.475,0.306,0.559]})
display(df)
print("  ")


# In[91]:



print (" ")
print ("Frequency = 200 MHz")

df = pd.DataFrame({"Tissue":["Air","Skin","Fat","Muscle",
                             "Cartilage","Lung","Brain(GM)"],
                             "Epsilon(Permittivity)":
                             [1,55.7,12.0,60.2,49.2,26.6,65.1],
                             "Sigma(Conductivity)":
                             [0.001,0.582,0.0726,0.743,0.518,0.335,0.639]})
display(df)
print("  ")


# In[92]:



print (" ")
print ("Frequency = 300 MHz")

df = pd.DataFrame({"Tissue":["Air","Skin","Fat","Muscle",
                             "Cartilage","Lung","Brain(GM)"],
                           "Epsilon(Permittivity)":
                           [1,49.8,11.7,58.2,46.8,24.8,60.0],
                           "Sigma(Conductivity)":
                           [0.001,0.641,0.0765,0.771,0.553,0.356,0.692]})
display(df)
print("  ")


# In[93]:



print (" ")
print ("Frequency = 400 MHz")

df = pd.DataFrame({"Tissue":["Air","Skin","Fat","Muscle",
                             "Cartilage","Lung","Brain(GM)"],
                             "Epsilon(Permittivity)":
                             [1,46.8,11.6,57.1,45.5,23.8,57.4],
                             "Sigma(Conductivity)":
                             [0.001,0.688,0.0807,0.796,0.586,0.374,0.737]})
display(df)
print("  ")


# In[94]:


print (" ")
print ("Frequency = 500 MHz")

df = pd.DataFrame({"Tissue":["Air","Skin","Fat","Muscle",
                             "Cartilage","Lung","Brain(GM)"],
                             "Epsilon(Permittivity)":
                             [1,44.9,11.5,56.4,44.6,23.2,55.8],
                             "Sigma(Conductivity)":
                             [0.001,0.728,0.0854,0.822,0.621,0.391,0.779]})
display(df)
print("  ")


# In[95]:



print (" ")
print ("Frequency = 600 MHz")

df = pd.DataFrame({"Tissue":["Air","Skin","Fat","Muscle",
                             "Cartilage","Lung","Brain(GM)"],
                             "Epsilon(Permittivity)":
                             [1,43.6,11.5,56.0,44.0,22.8,54.7],
                             "Sigma(Conductivity)":
                             [0.001,0.765,0.0905,0.850,0.658,0.407,0.819]})
display(df)
print("  ")


# In[96]:



print (" ")
print ("Frequency = 700 MHz")

df = pd.DataFrame({"Tissue":["Air","Skin","Fat","Muscle",
                             "Cartilage","Lung","Brain(GM)"],
                             "Epsilon(Permittivity)":
                             [1,42.7,11.4,55.6,43.5,22.5,53.9],
                             "Sigma(Conductivity)":
                             [0.001,0.800,0.0962,0.879,0.697,0.423,0.860]})
display(df)
print("  ")


# In[97]:



print (" ")
print ("Frequency = 800 MHz")

df = pd.DataFrame({"Tissue":["Air","Skin","Fat","Muscle",
                             "Cartilage","Lung","Brain(GM)"],
                             "Epsilon(Permittivity)":
                             [1,42.0,11.4,55.3,43.0,22.2,53.3],
                             "Sigma(Conductivity)":
                             [0.001,0.834,0.102,0.910,0.738,0.440,0.900]})

display(df)
print("  ")


# In[ ]:
