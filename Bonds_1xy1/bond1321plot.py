# The below code was created by Biswajit Sadhu

import math
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import font_manager
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)

# Coefficients for iso-pentane (in kcal/mol) (UA 1999 paper: J. Phys. Chem. #B1999,103,4508-4517
c0 = -0.4992
c1 = 0.8525
c2 = -0.2224
c3 = 0.8774

# Total_torsion is the function that produces the energy landscape
total_torsion = []
energy_x_list = []
energy_y_list = []
# Goes from 0 to 360 degrees
angles = np.arange(0,360,0.1)
for phi in np.arange(0,360,0.1):
    phi = np.deg2rad(phi)
    energy_x = c0 + c1*(1+np.cos(phi-np.deg2rad(56))) + c2*(1-np.cos(2*(phi-np.deg2rad(56))))+ c3*(1+np.cos(3*(phi-np.deg2rad(56))))
    energy_y = c0 + c1*(1+np.cos(phi+np.deg2rad(56))) + c2*(1-np.cos(2*(phi+np.deg2rad(56))))+ c3*(1+np.cos(3*(phi+np.deg2rad(56))))
    e = energy_x + energy_y
    total_torsion.append(e)
    energy_x_list.append(energy_x)
    energy_y_list.append(energy_y)

# Plots the function; all of the display specifications
fig, (ax1) = plt.subplots(figsize=(10, 6),ncols=1)
matplotlib.rcParams.update({'font.size': 16})
matplotlib.rcParams['font.family'] = ['DejaVu Serif']
font = font_manager.FontProperties(family='DejaVu Serif',
                                   weight='normal',
                                   style='normal', size=16)
ax1.set_ylim(-1,5,1)
ax1.set_xlim(0,360,60)
ax1.set_xticks((np.arange(0,360,60)))
ax1.set_yticks((np.arange(0,6,1)))
#ax1.axes.xaxis.set_visible(False)
#ax1.axes.yaxis.set_visible(False)
ax1.minorticks_on()
ax1.tick_params(which='minor', length=4, color='black')
ax1.plot(angles,total_torsion,color='orange',linewidth=3)
#ax1.plot(angles,energy_x_list,color='black',linewidth=3,linestyle='dashed')
#ax1.plot(angles,energy_y_list,color='black',linewidth=4,linestyle='dotted')
ax1.set_xlabel('Dihedral angle ($\degree$)', fontsize=20)
ax1.set_ylabel('Energy (kcal/mol)', fontsize=20)

fig.tight_layout()
ax1.legend(['CH$_3$-CH-CH-(CH$_3$)$_2$'], frameon=False)
plt.savefig('plot1321.png', dpi=600)
