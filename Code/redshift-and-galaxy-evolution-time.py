#!/usr/bin/env python
# coding: utf-8

# In[1]:


import math
import matplotlib.pyplot as plt
import numpy as np

Omega_L0 = 0.6911
Omega_r0 = 0
Omega_m0 = 0.3089
h_0 = 0.6774
#Note the difference between h_0 here and h below 
#define lookback time[Gyr] using the second method
n_integral = 3000              
def f(z):
    s = 1/((1+z)*math.sqrt(((1+z)**2)*(1+Omega_m0*z)-z*(2+z)*Omega_L0))
    return s
def lookbacktime(z):
    h = (z-0)/(2.0000000*n_integral)
    F0 = f(0)+f(z)
    F1 = 0
    F2 = 0
    for j in range(1,2*n_integral):
        x = j*h
        if j%2 == 0:
            F2 = F2+f(x)
        else:
            F1 = F1+f(x)
    return ((h/3)*(F0+(2*F2)+(4*F1)))*9.78*(10**9)/h_0/(10**9)#[Gyr]

#get an array of Z
Redshift_raw = dict(np.load('/huawei/osv1/chenyaoxin/workspace/Data/SnapRedshift.npz'))
#data from TNG300
z_values_raw = np.append(np.delete(np.append(Redshift_raw['Redshift'][13:99],[0])[::-1],[86]),[6.0])
t_values_raw = np.array([lookbacktime(6)-lookbacktime(z) for z in z_values_raw])
delta_t=10**(-5)

flag = 0
z_values = np.array([0])
while flag <= (z_values_raw.shape[0]-2):
    dim = z_values.shape[0]
    z_values = np.delete(z_values, [dim-1])
    n = int((abs(t_values_raw[flag]-t_values_raw[flag+1]))/delta_t)
    z_values = np.append(z_values, np.linspace(z_values_raw[flag],z_values_raw[flag+1],n))
    flag += 1

t_values = np.array([lookbacktime(6)-lookbacktime(z) for z in z_values])
np.savez('/huawei/osv1/chenyaoxin/workspace/Data/redshift-evolution_time10**(-5).npz', redshift = z_values, evolution_time = t_values)


# In[ ]:




