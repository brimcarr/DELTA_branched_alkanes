# The below code was created by Biswajit Sadhu

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import font_manager
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

#CH3-CH-CH-CH3(in kcal/mol) (1-3-3-1)
#The total number of dihedrals around each rotatable bond is
#given by the product of the number of attached groups on one end
#of the bond and the number of attached groups on the other end.

c0 = -0.4992
c1 = 0.8525
c2 = -0.2224
c3 = 0.8774

energy_x1_list = []
energy_x2_list = []
energy_x3_list = []
energy_x4_list = []
#sum of x1 and x2 (practically resembles to the case of 1-2-3-1)
energy_one_torsion = []
#sum of all
total_torsion = []

angles = np.arange(0,360,0.1)
for phi in np.arange(0,360,0.1):
    phi = np.deg2rad(phi)
    #CH3-CH-CH-CH3 TAKE THE BISECTOR AS THE DIHEDRAL AXIS
    mphi = (phi+np.deg2rad(-56))
    energy_x1 = [c0
                 + c1*(1+np.cos(1*mphi))
                 + c2*(1+np.cos((2*mphi)-np.deg2rad(180)))
                 + c3*(1+np.cos(3*mphi))]
    mphi = (phi+np.deg2rad(56))
    energy_x2 = [c0
                 + c1*(1+np.cos(1*mphi))
                 + c2*(1+np.cos((2*mphi)-np.deg2rad(180)))
                 + c3*(1+np.cos(3*mphi))]

    mphi = (phi+np.deg2rad(-56-112))
    energy_x3 = [c0
                 + c1*(1+np.cos(1*mphi))
                 + c2*(1+np.cos((2*mphi)-np.deg2rad(180)))
                 + c3*(1+np.cos(3*mphi))]

    mphi = (phi+np.deg2rad(56-112))
    energy_x4 = [c0
                 + c1*(1+np.cos(1*mphi))
                 + c2*(1+np.cos((2*mphi)-np.deg2rad(180)))
                 + c3*(1+np.cos(3*mphi))]

    e = energy_x1[0] + energy_x2[0] + energy_x3[0] + energy_x4[0]

    #e = energy_x1[0] + energy_y1[0] + energy_z1[0]
    total_torsion.append(e)
    energy_x1_list.append(energy_x1[0])
    energy_x2_list.append(energy_x2[0])
    energy_x3_list.append(energy_x3[0])
    energy_x4_list.append(energy_x4[0])
    energy_one_torsion.append(energy_x1[0]+energy_x2[0])

fig, (ax1) = plt.subplots(figsize=(10, 6),ncols=1)
matplotlib.rcParams.update({'font.size': 16})
matplotlib.rcParams['font.family'] = ['DejaVu Serif']
font = font_manager.FontProperties(family='DejaVu Serif',
                                   weight='normal',
                                   style='normal', size=16)
# Compute and print the critical values.
print('abs min (alpha)=',min(total_torsion),'local min (beta)=',total_torsion[1200],
'local max (gamma)=',total_torsion[1800],'abs max (delta)=',max(total_torsion))
#ax1.set_ylim(-10,20,3)
ax1.set_xlim(0,360,60)
ax1.set_xticks((np.arange(0,370,60)))
#ax1.set_yticks((np.arange(-10,20,3)))
ax1.minorticks_on()
ax1.tick_params(which='minor', length=4, color='black')
ax1.plot(angles,total_torsion,color='#ffd800',linewidth=3)
#ax1.plot(angles,energy_one_torsion,color='grey',linewidth=3,linestyle=None)
ax1.plot(angles,energy_x1_list,color='black',linewidth=3,linestyle='dashed')
ax1.plot(angles,energy_x2_list,color='black',linewidth=4,linestyle='dotted')
ax1.plot(angles,energy_x3_list,color='grey',linewidth=3,linestyle='dashed')
ax1.plot(angles,energy_x4_list,color='grey',linewidth=4,linestyle='dotted')
#ax1.hlines(0,0,360,color='grey',linewidth=1,linestyle='dotted')

#ax1.plot(angles,energy_z1_list,color='purple',linewidth=3,linestyle='-.')
ax1.set_xlabel('Dihedral angle ($\degree$)', fontsize=20)
ax1.set_ylabel('Energy (kcal/mol)', fontsize=20)


fig.tight_layout()
ax1.legend(['(CH$_3$)$_2$-CH-CH-(CH$_3$)$_2$'], frameon=False)
#ax1.legend(['(CH$_3$)$_2$-CH-CH-(CH$_3$)$_2$','CH$_3$-CH-CH-(CH$_3$)$_2$'], frameon=False)
plt.savefig('plot1331.png', dpi=600)
