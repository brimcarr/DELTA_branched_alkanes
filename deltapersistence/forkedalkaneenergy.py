"""Methods for computing and working with forked alkane energy landscapes.

The potential energy of the forked alkane molecules admits an approximate
analytical formula. This module computes that potential energy directly
at a given point or over a mesh. It also implements a more sophisticated
analytical formula for the cases of butane and pentane. Several other
utilities for working with forked alkane molecule names are provided.
Modified from code written by Joshua Mirth!


    Copyright (C) 2020 Joshua Mirth, modified by Brittany Story 2022
    (Added new functions to expand the use of Mirth's code)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.


Functions
---------

* alkane_to_carbons
* carbons_to_alkane
* n_alkane_energy_mesh
* nonbonded_potential
* potential

"""

import numpy as np

# Set constants and rename things for simplicity.
# Numpy operations know to act component-wise on arrays.
PI = np.pi
cos = np.cos
sin = np.sin
sqrt = np.sqrt

def n_alkane_energy_mesh(bondlist, n, res=.1, nonbonded=False):
    """Create a mesh of potential energy values.

    Meshes the cube [0,2pi]^n and computes values of the potential
    energy of the alkane on that mesh.

    Parameters
    ----------
    n : int
        Dimension of the array (n = number of carbons - 3).
    res : float, optional
        Resolution at which to mesh the grid. Default is 0.1, which is
        reasonable for most examples.
    nonbonded : bool, optional
        If true, use the more detailed formula which accounts for non-
        bonded potentials. Default is false. Currently true is only
        available if n=2 or n=3.

    Returns
    -------
    V : numpy ndarray
        Mesh of potential energy function values.

    """
    mesh = n*(np.arange(0,2*PI,res),)
    phi = np.array(np.meshgrid(*mesh, indexing='xy'))
    if nonbonded:
        V = nonbonded_potential(phi)
    else:
        V = potential_ua(phi,bondlist)
    return V


def potential_aa(phi):
    """Simple analytical formula for the OPLS-aa n-alkane potential energy.

    Formula derived by members of YZ's group. This does not account for
    the energy between non-bonded carbon atoms.

    Parameters
    ----------
    phi : numpy ndarray
        Mesh of the domain (n copies of [0,pi]).

    Returns
    -------
    V : numpy ndarray
        Function values over the domain.

    """

    V = ((2.9288-1.4644*cos(phi-PI)
        + 0.2092*cos(phi-PI)**2-1.6736*cos(phi-PI)**3)
        + 0.6276+1.8828*cos(phi-PI)-2.5104*cos(phi-PI)**3
        + 3*(0.6276+1.8828*cos(-phi-(1/3)*PI)-2.5104*cos(-phi-(1/3)*PI)**3)
        + 3*(0.6276+1.8828*cos(phi-(1/3)*PI)-2.5104*cos(phi-(1/3)*PI)**3))
    V = sum(V[:])      # V = V(phi_1) + V(phi_2) + ...
    return V


def potential_ua(phi,bondlist):
    #Brittany is editing this part.
    # Need list(n1,n2,n3 etc.)
    """Simple analytical formula for the OPLS-ua branched alkane potential energy.

    Formulas derived by members of Aurora's group. This does not account for
    the energy between non-bonded carbon atoms.

    Parameters
    ----------
    phi : numpy ndarray
        Mesh of the domain (n copies of [0,pi]).

    Returns
    -------
    V : numpy ndarray
        Function values over the domain.

    """
# This is the function for the energy landscape of dihedral type 1-2-2-1
    c10 = 0
    c11 = 0.7059
    c12 = -0.1355
    c13 = 1.5735
    n1 = bondlist[0]
    b1phi = phi[0:n1]

    V1 = (c10 + c11*(1+np.cos(b1phi))
                + c12*(1-np.cos(2*(b1phi)))
                + c13*(1+np.cos(3*(b1phi))))
    V1 = sum(V1[:])      # V1 = V1(phi_1) + V1(phi_2) + ... + V1(phi_{n_1-1})

# This is the function for the energy landscape of dihedral type 1-3-2-1 (isopentane)
    c20 = -0.4992
    c21 = 0.8525
    c22 = -0.2224
    c23 = 0.8774
    n2 = bondlist[1]
    b2phi = phi[n1:n1+n2]

    energy_x21 = (c20 + c21*(1+np.cos(b2phi-np.deg2rad(56)))
                + c22*(1-np.cos(2*(b2phi-np.deg2rad(56))))
                + c23*(1+np.cos(3*(b2phi-np.deg2rad(56)))))
    energy_x22 = (c20 + c21*(1+np.cos(b2phi+np.deg2rad(56)))
                + c22*(1-np.cos(2*(b2phi+np.deg2rad(56))))
                + c23*(1+np.cos(3*(b2phi+np.deg2rad(56)))))
    V2 = energy_x21 + energy_x22
    V2 = sum(V2[:])      # V = V(phi_1) + V(phi_2) + ...

# This is the function for the energy landscape of dihedral type 1-3-3-1
    c30 = -0.4992
    c31 = 0.8525
    c32 = -0.2224
    c33 = 0.8774
    n3 = bondlist[2]
    b3phi = phi[n1+n2:n1+n2+n3]

    mphi = (b3phi+np.deg2rad(-56))
    energy_x31 = [c30
                 + c31*(1+np.cos(1*mphi))
                 + c32*(1+np.cos((2*mphi)-np.deg2rad(180)))
                 + c33*(1+np.cos(3*mphi))]
    mphi = (b3phi+np.deg2rad(56))
    energy_x32 = [c30
                 + c31*(1+np.cos(1*mphi))
                 + c32*(1+np.cos((2*mphi)-np.deg2rad(180)))
                 + c33*(1+np.cos(3*mphi))]

    mphi = (b3phi+np.deg2rad(-56-112))
    energy_x33 = [c30
                 + c31*(1+np.cos(1*mphi))
                 + c32*(1+np.cos((2*mphi)-np.deg2rad(180)))
                 + c33*(1+np.cos(3*mphi))]

    mphi = (b3phi+np.deg2rad(56-112))
    energy_x34 = [c30
                 + c31*(1+np.cos(1*mphi))
                 + c32*(1+np.cos((2*mphi)-np.deg2rad(180)))
                 + c33*(1+np.cos(3*mphi))]

    V3 = energy_x31[0] + energy_x32[0] + energy_x33[0] + energy_x34[0]
    V3 = sum(V3[:])      # V = V(phi_1) + V(phi_2) + ...

# This is the function for the energy landscape of dihedral type 1-4-2-1
    c40 = 0
    c41 = 0
    c42 = 0
    c43 = 0.9172

    n4 = bondlist[3]
    b4phi = phi[n1+n2+n3:n1+n2+n3+n4]

    mphi = (b4phi+np.deg2rad(-(107/2) - 107))

    energy_x40 = [c40
                 + c41*(1+np.cos(1*mphi))
                 + c42*(1+np.cos((2*mphi)-np.deg2rad(180)))
                 + c43*(1+np.cos(3*mphi))]

    mphi = (b4phi+np.deg2rad(-107/2))

    energy_x41 = [c40
                 + c41*(1+np.cos(1*mphi))
                 + c42*(1+np.cos((2*mphi)-np.deg2rad(180)))
                 + c43*(1+np.cos(3*mphi))]

    mphi = (b4phi+np.deg2rad(107/2))

    energy_x42 = [c40
                 + c41*(1+np.cos(1*mphi))
                 + c42*(1+np.cos((2*mphi)-np.deg2rad(180)))
                 + c43*(1+np.cos(3*mphi))]

    mphi = (b4phi+np.deg2rad(-(107/2)-107-107))

    V4 = energy_x40[0] + energy_x41[0] + energy_x42[0]

    V4 = sum(V4[:])      # V = V(phi_1) + V(phi_2) + ...

# This is the function for the energy landscape of dihedral type 1-4-3-1
    c50 = 0
    c51 = 0
    c52 = 0
    c53 = 2.7221
    n5 = bondlist[4]
    b5phi = phi[n1+n2+n3+n4:n1+n2+n3+n4+n5]

    #CH3-C-C-CH3 TAKE THE BISECTOR of two methyl group AS THE DIHEDRAL AXIS
    #for the methyl with opposite to dihedral axis
    mphi = (b5phi+np.deg2rad(-(107/2) - 107))

    energy_x50 = [c50
                + c51*(1+np.cos(1*mphi))
                + c52*(1+np.cos((2*mphi)-np.deg2rad(180)))
                + c53*(1+np.cos(3*mphi))]

    mphi = (b5phi+np.deg2rad(-107/2))

    energy_x51 = [c50
                + c51*(1+np.cos(1*mphi))
                + c52*(1+np.cos((2*mphi)-np.deg2rad(180)))
                + c53*(1+np.cos(3*mphi))]

    mphi = (b5phi+np.deg2rad(107/2))

    energy_x52 = [c50
                + c51*(1+np.cos(1*mphi))
                + c52*(1+np.cos((2*mphi)-np.deg2rad(180)))
                + c53*(1+np.cos(3*mphi))]

    mphi = (b5phi+np.deg2rad(-(107/2)-107-107))

    energy_x53 = [c50
                + c51*(1+np.cos(1*mphi))
                + c52*(1+np.cos((2*mphi)-np.deg2rad(180)))
                + c53*(1+np.cos(3*mphi))]

    mphi = (b5phi+np.deg2rad(-(107/2)-107))

    energy_x54 = [c50
                + c51*(1+np.cos(1*mphi))
                + c52*(1+np.cos((2*mphi)-np.deg2rad(180)))
                + c53*(1+np.cos(3*mphi))]

    mphi = (b5phi+np.deg2rad((107/2)-107))

    energy_x55 = [c50
                + c51*(1+np.cos(1*mphi))
                + c52*(1+np.cos((2*mphi)-np.deg2rad(180)))
                + c53*(1+np.cos(3*mphi))]

    V5 = energy_x50[0] + energy_x51[0] + energy_x52[0] + energy_x53[0] + energy_x54[0] + energy_x55[0]
    V5 = sum(V5[:])      # V = V(phi_1) + V(phi_2) + ...

# This is the function for the energy landscape of dihedral type 1-4-4-1
    c60 = 5.8677
    c61 = 0
    c62 = 0
    c63 = -4.0704
    c66 = 1.0523
    c612 = 0.2985
    n6 = bondlist[5]
    b6phi = phi[n1+n2+n3+n4+n5:n1+n2+n3+n4+n5+n6]
    mphi = (b6phi+np.deg2rad(-(107/2) - 107))

    energy_x60 = [c60
                 + c63*(1+np.cos(3*mphi))
                 + c66*(1+np.cos((6*mphi)-np.deg2rad(180)))
                 + c612*(1+np.cos(12*mphi))]

    mphi = (b6phi+np.deg2rad(-107/2))

    energy_x61 = [c60
                 + c63*(1+np.cos(3*mphi))
                 + c66*(1+np.cos((6*mphi)-np.deg2rad(180)))
                 + c612*(1+np.cos(12*mphi))]

    mphi = (b6phi+np.deg2rad(107/2))

    energy_x62 = [c60
                 + c63*(1+np.cos(3*mphi))
                 + c66*(1+np.cos((6*mphi)-np.deg2rad(180)))
                 + c612*(1+np.cos(12*mphi))]

    mphi = (b6phi+np.deg2rad(-(107/2)-107-107))

    energy_x63 = [c60
                 + c63*(1+np.cos(3*mphi))
                 + c66*(1+np.cos((6*mphi)-np.deg2rad(180)))
                 + c612*(1+np.cos(12*mphi))]

    mphi = (b6phi+np.deg2rad(-(107/2)-107))

    energy_x64 = [c60
                 + c63*(1+np.cos(3*mphi))
                 + c66*(1+np.cos((6*mphi)-np.deg2rad(180)))
                 + c612*(1+np.cos(12*mphi))]

    mphi = (b6phi+np.deg2rad((107/2)-107))

    energy_x65 = [c60
                 + c63*(1+np.cos(3*mphi))
                 + c66*(1+np.cos((6*mphi)-np.deg2rad(180)))
                 + c612*(1+np.cos(12*mphi))]

    mphi = (b6phi+np.deg2rad(-(107/2)-107+107))

    energy_x66 = [c60
                 + c63*(1+np.cos(3*mphi))
                 + c66*(1+np.cos((6*mphi)-np.deg2rad(180)))
                 + c612*(1+np.cos(12*mphi))]

    mphi = (b6phi+np.deg2rad(-(107/2)+107))

    energy_x67 = [c60
                 + c63*(1+np.cos(3*mphi))
                 + c66*(1+np.cos((6*mphi)-np.deg2rad(180)))
                 + c612*(1+np.cos(12*mphi))]

    mphi = (b6phi+np.deg2rad((107/2)+107))

    energy_x68 = [c60
                 + c63*(1+np.cos(3*mphi))
                 + c66*(1+np.cos((6*mphi)-np.deg2rad(180)))
                 + c612*(1+np.cos(12*mphi))]

    V6 = energy_x60[0] + energy_x61[0] + energy_x62[0] + energy_x63[0] + energy_x64[0]
    + energy_x65[0] + energy_x66[0] + energy_x67[0] + + energy_x68[0]

    V6 = sum(V6[:])      # V = V(phi_1) + V(phi_2) + ...
    return V1 + V2 + V3 + V4 + V5 + V6

def nonbonded_potential(phi):
    """Detailed analytical formula for n-alkane potential energy.

    Formula derived by members of YZ's group. This version accounts for
    interactions between non-bonded atoms. Currently only possible for
    dimensions 2 and 3 (pentane and hexane).

    Parameters
    ----------
    phi : numpy ndarray
        Mesh of the domain (n copies of [0,pi]).

    Returns
    -------
    V : numpy ndarray
        Function values over the domain.

    Raises
    ------
    ValueError
        If the dimension of the input domain array is not 2 or 3.

    """

    # Physical constants.
    c1 = 2.9519
    c2 = -0.5670
    c3 = 6.5794
    r0 = 0.1529
    eps = 0.7322
    sig = 0.3905
    n = len(phi)
    if n == 2:  # Pentane formula.
        r15 = ((4*sqrt(3)/9)*r0*sqrt(11-4*(cos(phi[0])+cos(phi[1]))
                - 2*cos(phi[0]+phi[1])+cos(phi[0]-phi[1])))
        Vd = (c1*(1+cos(phi[0]))+c2*(1-cos(2*phi[0]))+c3*(1+cos(3*phi[0]))
              + c1*(1+cos(phi[1]))+c2*(1-cos(2*phi[1]))+c3*(1+cos(3*phi[1])))
        Vnb = (4*eps*((sig/r15)**12 - (sig/r15)**6))
        V = Vd + Vnb
        return V
    elif n == 3:    # Hexane formula.
        Vd = (c1*(1+cos(phi[0]))+c2*(1-cos(2*phi[0]))+c3*(1+cos(3*phi[0]))
              + c1*(1+cos(phi[1]))+c2*(1-cos(2*phi[1]))+c3*(1+cos(3*phi[1]))
              + c1*(1+cos(phi[2]))+c2*(1-cos(2*phi[2]))+c3*(1+cos(3*phi[2])))
        d21 = (sqrt((48/81)*r0*r0*(11-4*(cos(phi[0])+cos(phi[1]))
               - 2*cos(phi[0]+phi[1])+cos(phi[0]-phi[1]))))
        d22 = (sqrt((48/81)*r0*r0*(11-4*(cos(phi[1])+cos(phi[2]))
               - 2*cos(phi[1]+phi[2])+cos(phi[1]-phi[2]))))
        d3 = (sqrt((1/81)*r0*r0*(689-208*(cos(phi[0])+cos(phi[2]))
               - 256*cos(phi[1])+192*sin(phi[1])*(sin(phi[0])+sin(phi[2]))
               - 64*cos(phi[1])*(cos(phi[0])+cos(phi[2]))
               + 128*cos(phi[0])*cos(phi[2])
               + 48*sin(phi[1])*(cos(phi[0])*sin(phi[2])
               - sin(phi[0])*cos(phi[2]))
               + 16*cos(phi[1])*(9*sin(phi[0])*sin(phi[2])
               - cos(phi[0])*cos(phi[2])))))
        Vnb = (4*eps*((sig/d21)**12-(sig/d21)**6)+4*eps*((sig/d22)**12
               - (sig/d22)**6)+4*eps*((sig/d3)**12-(sig/d3)**6))
        V = Vd + Vnb;
        return V
    else:
        raise ValueError('For non-bonded energy, n must be 2 or 3 (pentane or '
                         'hexane).\n The value of n was: {}'.format(n))

def alkane_to_carbons(alkane):
    """Convert an alkane name to the number of carbons it contains.

    Parameters
    ----------
    alkane : string
        Name of an n-alkane molecule

    Returns
    -------
    carbon : int
        Number of carbons in the n-alkane

    """

    carbon_dict = {
        'methane': 1,
        'ethane': 2,
        'propane': 3,
        'butane': 4,
        'pentane': 5,
        'hexane': 6,
        'heptane': 7,
        'octane': 8,
        'nonane': 9,
        'decane': 10,
        'undecane': 11,
        'dodecane': 12,
        'tridecane': 13,
        'tetradecane': 14,
        'pentadecane': 15,
        'hexadecane': 16,
        'heptadecane': 17,
        'octadecane': 18,
        'nonadecane': 19,
        'icosane': 20,
        'henicosane': 21,
        'docosane': 22,
        'tricosane': 23,
        'tetracosane': 24,
        'pentacosane': 25,
        'hexacosane': 26,
        'heptacosane': 27,
        'octacosane': 28,
        'nonacosane': 29,
        'triacontane': 30,
        'hentriacontane': 31,
        'dotriacontane': 32,
        'tritriacontane': 33,
        'tetratriacontane': 34,
        'pentatriacontane': 35,
        'hexatriacontane': 36,
        'heptatriacontane': 37,
        'octatriacontane': 38,
        'nonatriacontane': 39,
        'tetracontane': 40,
        'hentetracontane': 41,
        'dotetracontane': 42,
        'tritetracontane': 43,
        'tetratetracontane': 44,
        'pentatetracontane': 45,
        'hexatetracontane': 46,
        'heptatetracontane': 47,
        'octatetracontane': 48,
        'nonatetracontane': 49,
        'pentacontane': 50,
        'henpentacontane': 51,
        'dopentacontane': 52,
        'tripentacontane': 53,
        'tetrapentacontane': 54,
        'pentapentacontane': 55,
        'hexapentacontane': 56,
        'heptapentacontane': 57,
        'octapentacontane': 58,
        'nonapentacontane': 59,
        'hexacontane': 60,
        'henhexacontane': 61,
        'dohexacontane': 62,
        'trihexacontane': 63,
        'tetrahexacontane': 64,
        'pentahexacontane': 65,
        'hexahexacontane': 66,
        'heptahexacontane': 67,
        'octahexacontane': 68,
        'nonahexacontane': 69,
        'heptacontane': 70,
        'henheptacontane': 71,
        'doheptacontane': 72,
        'triheptacontane': 73,
        'tetraheptacontane': 74,
        'pentaheptacontane': 75,
        'hexaheptacontane': 76,
        'heptaheptacontane': 77,
        'octaheptacontane': 78,
        'nonaheptacontane': 79,
        'octacontane': 80,
        'henoctacontane': 81,
        'dooctacontane': 82,
        'trioctaane': 83,
        'tetraoctacontane': 84,
        'pentaoctacontane': 85,
        'hexaoctacontane': 86,
        'heptaoctacontane': 87,
        'octaoctacontane': 88,
        'nonaoctacontane': 89,
        'nonacontane': 90,
        'hennonacontane': 91,
        'dononacontane': 92,
        'trinonacontane': 93,
        'tetranonacontane': 94,
        'pentanonacontane': 95,
        'hexanonacontane': 96,
        'heptanonacontane': 97,
        'octanonacontane': 98,
        'nonanonacontane': 99,
        'hectane': 100,
        'henihectane': 101,
        'dohectane': 102,
        'trihectane': 103,
        'tetrahectane': 104,
        'pentahectane': 105,
        'hexahectane': 106,
        'heptahectane': 107,
        'octahectane': 108,
        'nonahectane': 109,
        'decahectane': 110,
        'undecahectane': 111,
        'dodecahectane': 112,
        'tridecahectane': 113,
        'tetradecahectane': 114,
        'pentadecahectane': 115,
        'hexadecahectane': 116,
        'heptadecahectane': 117,
        'octadecahectane': 118,
        'nonadecahectane': 119,
        'icosahectane': 120
    }
    alkane = alkane.casefold()
    try:
        carbons = carbon_dict[alkane]
        return carbons
    except KeyError:
        print('The n-alkane molecule ' + alkane + ' is not recognized.')
        raise

def carbons_to_alkane(carbons):
    """Gives the name of the n-alkane molecule with a given number of
    carbon atoms.

    Parameters
    ----------
    carbons : int
        Number of carbons in the n-alkane.

    Returns
    -------
    alkane : string
        Name of the n-alkane molecule with `carbons` carbon atoms.

    """

    molecule_dict = {
        1: 'methane',
        2: 'ethane',
        3: 'propane',
        4: 'butane',
        5: 'pentane',
        6: 'hexane',
        7: 'heptane',
        8: 'octane',
        9: 'nonane',
        10: 'decane',
        11: 'undecane',
        12: 'dodecane',
        13: 'tridecane',
        14: 'tetradecane',
        15: 'pentadecane',
        16: 'hexadecan',
        17: 'heptadecane',
        18: 'octadecane',
        19: 'nonadecane',
        20: 'icosane',
        21: 'henicosane',
        22: 'docosane',
        23: 'tricosane',
        24: 'tetracosane',
        25: 'pentacosane',
        26: 'hexacosane',
        27: 'heptacosane',
        28: 'octacosane',
        29: 'nonacosane',
        30: 'triacontane',
        31: 'hentriacontane',
        32: 'dotriacontane',
        33: 'tritriacontane',
        34: 'tetratriacontane',
        35: 'pentatriacontane',
        36: 'hexatriacontane',
        37: 'heptatriacontane',
        38: 'octatriacontane',
        39: 'nonatriacontane',
        40: 'tetracontane',
        41: 'hentetracontane',
        42: 'dotetracontane',
        43: 'tritetracontane',
        44: 'tetratetracontane',
        45: 'pentatetracontane',
        46: 'hexatetracontane',
        47: 'heptatetracontane',
        48: 'octatetracontane',
        49: 'nonatetracontane',
        50: 'pentacontane',
        51: 'henpentacontane',
        52: 'dopentacontane',
        53: 'tripentacontane',
        54: 'tetrapentacontane',
        55: 'pentapentacontane',
        56: 'hexapentacontane',
        57: 'heptapentacontane',
        58: 'octapentacontane',
        59: 'nonapentacontane',
        60: 'hexacontane',
        61: 'henhexacontane',
        62: 'dohexacontane',
        63: 'trihexacontane',
        64: 'tetrahexacontane',
        65: 'pentahexacontane',
        66: 'hexahexacontane',
        67: 'heptahexacontane',
        68: 'octahexacontane',
        69: 'nonahexacontane',
        70: 'heptacontane',
        71: 'henheptacontane',
        72: 'doheptacontane',
        73: 'triheptacontane',
        74: 'tetraheptacontane',
        75: 'pentaheptacontane',
        76: 'hexaheptacontane',
        77: 'heptaheptacontane',
        78: 'octaheptacontane',
        79: 'nonaheptacontane',
        80: 'octacontane',
        81: 'henoctacontane',
        82: 'dooctacontane',
        83: 'trioctaane',
        84: 'tetraoctacontane',
        85: 'pentaoctacontane',
        86: 'hexaoctacontane',
        87: 'heptaoctacontane',
        88: 'octaoctacontane',
        89: 'nonaoctacontane',
        90: 'nonacontane',
        91: 'hennonacontane',
        92: 'dononacontane',
        93: 'trinonacontane',
        94: 'tetranonacontane',
        95: 'pentanonacontane',
        96: 'hexanonacontane',
        97: 'heptanonacontane',
        98: 'octanonacontane',
        99: 'nonanonacontane',
        100: 'hectane',
        101: 'henihectane',
        102: 'dohectane',
        103: 'trihectane',
        104: 'tetrahectane',
        105: 'pentahectane',
        106: 'hexahectane',
        107: 'heptahectane',
        108: 'octahectane',
        109: 'nonahectane',
        110: 'decahectane',
        111: 'undecahectane',
        112: 'dodecahectane',
        113: 'tridecahectane',
        114: 'tetradecahectane',
        115: 'pentadecahectane',
        116: 'hexadecahectane',
        117: 'heptadecahectane',
        118: 'octadecahectane',
        119: 'nonadecahectane',
        120: 'icosahectane'
    }
    try:
        molecule = molecule_dict[carbons]
        return molecule
    except KeyError:
        print('The name of the n-alkane molecule with ' + str(carbons) +
            ' carbon atoms is unknown. Must be an integer between 1 and 120.')
        raise
