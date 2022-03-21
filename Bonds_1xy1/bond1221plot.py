# The below code was created by Biswajit Sadhu

import math
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import font_manager
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)

# Coefficients for butane (in kcal/mol) (UA 1999 paper: J. Phys. Chem. #B1999,103,4508-4517
c0 = 0
c1 = 0.7059
c2 = -0.1355
c3 = 1.5735

total_torsion = []
# Goes from 0 to 360 degrees
angles = np.arange(0,360,0.1)
for phi in np.arange(0,360,0.1):
    phi = np.deg2rad(phi)
    e = c0 + c1*(1+np.cos(phi))+ c2*(1-np.cos(2*(phi)))+ c3*(1+np.cos(3*(phi)))

    total_torsion.append(e)

# print('abs min (alpha)=',min(total_torsion),'local min (beta)=',min(total_torsion[0:600]),
#        'local max (gamma)=',max(total_torsion[1800:2400]),'abs max (delta)=',max(total_torsion))


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
ax1.set_yticks((np.arange(0,5,1)))
ax1.axes.xaxis.set_visible(False)
ax1.axes.yaxis.set_visible(False)
ax1.minorticks_on()
ax1.tick_params(which='minor', length=4, color='black')
ax1.plot(angles,total_torsion,color='red',linewidth=3)
#ax1.set_xlabel('Dihedral angle ($\degree$)', fontsize=20)
#ax1.set_ylabel('Energy (kcal/mol)', fontsize=20)

fig.tight_layout()
#ax1.legend(['CH$_3$-CH-CH-CH$_3$'], frameon=False)
plt.savefig('plot1221.png', dpi=600)
