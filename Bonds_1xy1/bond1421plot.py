# The below code was created by Biswajit Sadhu

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import font_manager
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)


#CH3-CH2-C-CH3
c0 = 0
c1 = 0
c2 = 0
c3 = 0.9172



energy_x0_list = []
energy_x1_list = []
energy_x2_list = []

total_torsion = []

angles = np.arange(0,360,0.1)
for phi in np.arange(0,360,0.1):
    phi = np.deg2rad(phi)
    #CH3-C-C-CH3 TAKE THE BISECTOR of two methyl group AS THE DIHEDRAL AXIS
    #for the methyl with opposite to dihedral axis
    mphi = (phi+np.deg2rad(-(107/2) - 107))

    energy_x0 = [c0
                 + c1*(1+np.cos(1*mphi))
                 + c2*(1+np.cos((2*mphi)-np.deg2rad(180)))
                 + c3*(1+np.cos(3*mphi))]

    mphi = (phi+np.deg2rad(-107/2))

    energy_x1 = [c0
                 + c1*(1+np.cos(1*mphi))
                 + c2*(1+np.cos((2*mphi)-np.deg2rad(180)))
                 + c3*(1+np.cos(3*mphi))]

    mphi = (phi+np.deg2rad(107/2))

    energy_x2 = [c0
                 + c1*(1+np.cos(1*mphi))
                 + c2*(1+np.cos((2*mphi)-np.deg2rad(180)))
                 + c3*(1+np.cos(3*mphi))]

    mphi = (phi+np.deg2rad(-(107/2)-107-107))

    e = energy_x0[0] + energy_x1[0] + energy_x2[0]


    total_torsion.append(e)
    energy_x0_list.append(energy_x0[0])
    energy_x1_list.append(energy_x1[0])
    energy_x2_list.append(energy_x2[0])

print('abs min (alpha)=',min(total_torsion),'local min (beta)=',total_torsion[1200],
'local max (gamma)=',total_torsion[1800],'abs max (delta)=',max(total_torsion))


fig, (ax1) = plt.subplots(figsize=(10, 6),ncols=1)
matplotlib.rcParams.update({'font.size': 16})
matplotlib.rcParams['font.family'] = ['DejaVu Serif']
font = font_manager.FontProperties(family='DejaVu Serif',
                                   weight='normal',
                                   style='normal', size=16)
ax1.set_ylim(-1,6,3)
ax1.set_xlim(0,360,60)
ax1.set_xticks((np.arange(0,370,60)))
#ax1.set_yticks((np.arange(-10,20,3)))
ax1.minorticks_on()
ax1.tick_params(which='minor', length=4, color='black')
ax1.plot(angles,total_torsion,color='green',linewidth=3)
ax1.plot(angles,energy_x0_list,color='black',linewidth=3,linestyle='-.')
ax1.plot(angles,energy_x1_list,color='black',linewidth=3,linestyle='dashed')
ax1.plot(angles,energy_x2_list,color='black',linewidth=4,linestyle='dotted')

ax1.set_xlabel('Dihedral angle ($\degree$)', fontsize=20)
ax1.set_ylabel('Energy (kcal/mol)', fontsize=20)


fig.tight_layout()
ax1.legend(['(CH$_3$)-CH$_2$-C-(CH$_3$)$_3$'], frameon=False)
plt.savefig('plot1421.png', dpi=600)
