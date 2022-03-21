"""Script for computing sublevelset persistence.

    Copyright (C) 2020 Joshua Mirth and Johnathan Bush,
        modified by Brittany Story 2022 to expand the types of molecules
        the code can analyze.

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

"""

# External libraries
import os
import argparse
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plot
import matplotlib.patches as mpatches
import gudhi
from gudhi.persistence_graphical_tools import __min_birth_max_death as mb_md
# GUDHI plot routines try to use LaTeX for titles. Set this to True if
# you want that behavior.
gudhi.persistence_graphical_tools._gudhi_matplotlib_use_tex=False

import deltapersistence.forkedalkaneenergy as ae
import deltapersistence.forkedsublevelpersistence as sp
import deltapersistence.forkeddatautils as du

parser = argparse.ArgumentParser(
    description='Script for computing sublevelset persistence. Uses the '\
        'GUDHI library to compute the sublevelset persistence of a given '\
        'energy landscape. The energy landscape can be given as an energy '\
        'function (specified by either a random sample of function values or '\
        'an evenly-spaced mesh of function values) or an n-alkane molecule '\
        'can be specified and the appropriate energy landscape computed on '\
        'the fly. Minimum required arguments are an input file and filetype, '\
        'or a molecule specification.')
parser.add_argument('-b', '--bars', action='store_true',
    help='If specified, display the barcodes.')
parser.add_argument('-c', '--field', type=int, default=11,
    help='Coefficient field for computing homology. Must be a prime.')
parser.add_argument('-d', '--dimension', type=int, default = -1,
    help='Homological dimension to return. Default is -1 (all dimensions).')
parser.add_argument('-dg', '--diagram', action='store_true',
    help='If specified, display the persistence diagram.')
parser.add_argument('-ds', '--disp', action='store_true',
    help='If specified, print the persistence intervals to standard output.')
parser.add_argument('-e', '--energy', type=str, default='',
    help='Path to input file containing function values for sampled data.')
parser.add_argument('-f', '--file', type=str, default='',
    help='Path to input file.')
# parser.add_argument('-m', '--molecule', type=str, default='',
#     help='If provided, generates the energy landscape of the given n-alkane '\
#         'molecule and uses that as input.')
parser.add_argument('-n', '--carbons', type=int, default=0,
    help='If provided, generates the energy landscape of the alkane molecule '\
        'with the given number of carbon atoms and uses that as input.')
parser.add_argument('-nbonds', '--bonds', nargs='+', type=int, default=[0,0,0,0,0,0],
    help='If provided, generates the energy landscape of the molecule '\
        'with the given number of 1221, 1321, 1331, 1421, 1431, and 1441'\
        'dihedral types and uses that as input.')
parser.add_argument('-nb', '--nonbonded', action='store_true',
    help='If generating an energy landscape, use this flag to specify that '\
        'the formula accounting for non-bonded interactions should be used. '\
        'This is only valid for butane or pentane currently.')
parser.add_argument('-of', '--output_file', type=str, default='intervals',
    help='File path to write the output intervals to. The same file name is '\
        ' used for a text list of intervals, a barcode plot, or a diagram. '\
        'The correct file extension is automatically appended.')
parser.add_argument('-pd', '--periodic_dims', type=str, default=[],
    help='If a mesh is given, this lists the dimensions to treat as '\
        'periodic. This should be a string of Booleans, e.g. 101. If '\
        'provided for alkane molecule input, use 0 to set all dimensions to '\
        'nonperiodic. Any other input results in all periodic dimensions. If '\
        'provided for sampled data, any nonzero input sets all boundaries to '\
        'periodic. For sampled data this option is only available in '\
        'dimension two and three.')
parser.add_argument('-r', '--resolution', type=float, default=0.1,
    help='If generating an energy landscape, use this as the resolution.')
parser.add_argument('-t', '--input_type', type=str, default='',
    help='Type of input file. Options are "sample" or "mesh".')
parser.add_argument('-T', '--title', type=str, default='',
    help='Title for plots.')
parser.add_argument('-w', '--write', action='store_true',
    help='If specified, write the persistence intervals to a file. Passing '\
        'this argument saves all specified output formats. Intervals are '\
        'saved as a `.pers` file, while barcodes and diagrams are saved as '\
        '`.pdf` files.')

args = parser.parse_args()

# Format and validate input
input_type = args.input_type.lower()

print('######################################################################')

# Default to processing a file if one is specified, regardless of other
# arguments.
if args.file != '':
    if not os.path.isfile(args.file):
        print('The input file ' + args.file + ' was not found.')
        exit()
    if args.title == '':
        title = 'Sublevelset Persistence for ' + args.file
    else:
        title = args.title
    if input_type == 'sample':
        print('...Generating Delaunay triangulation from file ' + args.file)
        if args.periodic_dims != []:
            periodic = True
        else:
            periodic = False
        # TODO: check validity of data.
        coords = np.loadtxt(args.file)
        energy = np.loadtxt(args.energy)
        lower_star_cplx = sp.lower_star_delaunay(coords,energy,periodic)
        print('...The resulting complex contains ' +
            str(lower_star_cplx.num_simplices()) + ' simplices.')
        print('...Computing persistent homology.')
        starttime = datetime.now()
        intervals = lower_star_cplx.persistence(
            homology_coeff_field=args.field,persistence_dim_max=True)
        print('...Persistence computation complete.')
        print('...Elapsed time: ' + str(datetime.now() - starttime))
    else:   # The other parameter option is 'mesh', which will be treated as
            # the default case.
        print('...Reading energy landscape from file ' + args.file)
        # Energy landscape meshes should be saved as numpy .npy files.
        try:
            fnx = np.load(args.file)
        except ValueError:
            print('The file found at ' + args.file + ' is not a valid ' +
                'Numpy array stored in .npy format.')
            exit()
        except FileNotFoundError:
            print('The file ' + args.file + ' could not be found.')
            exit()
        if fnx.ndim != len(args.periodic_dims):
            print('...Supplied periodic dimensions not valid. Proceeding ' +
                'assuming no dimensions are periodic.')
            periodic_dims = [False]*fnx.ndim
        else:
            periodic_dims = []
            for s in args.periodic_dims:
                if s == '0' or s == 'f':
                    periodic_dims.append(False)
                else:
                    periodic_dims.append(True)
        print('...Generating complex and computing persistent homology.')
        starttime = datetime.now()
        intervals = sp.mesh_sublevel_persistence(
            fnx,
            periodic_dims,
            args.field)
        print('...Persistence computation complete.')
        print('...Elapsed time: ' + str(datetime.now() - starttime))
#elif args.carbons > 0 or args.molecule != '' or np.count_nonzero(args.bonds) > 0:
elif args.carbons > 0 or np.count_nonzero(args.bonds) > 0:
    if args.carbons > 3:
        # molecule = ae.carbons_to_alkane(args.carbons)
        carbons = args.carbons
        nlist = [carbons, 0, 0, 0, 0, 0]
    # elif args.molecule != '':
    #     molecule = args.molecule
    #     carbons = ae.alkane_to_carbons(args.molecule)
    elif np.count_nonzero(args.bonds) > 0:
        b1221 = args.bonds[0]
        b1321 = args.bonds[1]
        b1331 = args.bonds[2]
        b1421 = args.bonds[3]
        b1431 = args.bonds[4]
        b1441 = args.bonds[5]
        nlist = [b1221, b1321, b1331, b1421, b1431, b1441]
        carbons = sum(nlist) + 3
    else:
        print('Cannot compute the energy landscape of the branched alkane ' +
            'molecule with ' + str(args.carbons) + ' carbon atoms. ' +
            'The number of carbon atoms must be between 4 and 120.')
        exit()
    n = carbons - 3
    # Dimension of the energy landscape = number of bonds.
    # Dimension of the energy landscape = number of carbons - 3.
    bondlist = nlist

    if args.title == '':
    #     title = molecule.capitalize()
    # else:
        title = args.title
    print('...Computing the sublevelset persistence of ', n, 'carbon branched alkane.')
    if args.nonbonded:
        nb = 'taking'
    else:
        nb = 'not taking'
    print('...Determining the energy landscape of ' , n,  'carbon branched alkane to ' +
        'a resolution of ' + str(args.resolution) + ' and ' + nb + ' into ' +
        'account interactions between non-bonded carbons.')
    fnx_elements = int((2*np.pi/args.resolution)**n)
    print('...The function mesh has approximately ' + str(fnx_elements) +
        ' elements in it. If this is larger than 10^6 the ' +
        'computation may take a long time.')
    energy_landscape = ae.n_alkane_energy_mesh(
        bondlist,
        n,
        args.resolution,
        args.nonbonded)
    if args.periodic_dims == '0':
        periodic_dims = [False]*n
    else:
        periodic_dims = [True]*n
    print('...Computing persistent homology.')
    starttime = datetime.now()
    intervals = sp.mesh_sublevel_persistence(
        energy_landscape,
        periodic_dims,
        args.field)
    print('...Persistence computation complete.')
    print('...Elapsed time: ' + str(datetime.now() - starttime))

else:
    print('No valid input parameters provided. Run this script with the ' +
        'flag "--help" to see available input options.')
    exit()

if args.dimension >= 0:
    intervals = du.intervals_in_dimension(intervals,args.dimension)

if args.disp:
    print('Persistence intervals:')
    for i in intervals:
        print(i)
    if  args.write:
        interval_file = args.output_file + '.pers'
        print('Writing intervals to ' + interval_file)
        du.write_intervals(intervals, interval_file)

if args.bars:
    print('...Generating barcodes.')
    clean_intervals = du.remove_inf_bars(intervals)
    barcode_plot = gudhi.plot_persistence_barcode(
        clean_intervals,
        max_intervals=0,
        alpha=1.0)
    # Set formatting choices for barcodes.
    barcode_plot.set_title(title, fontsize='30', fontfamily='serif')
    barcode_plot.tick_params(labelsize='25')
    for tick in barcode_plot.get_xticklabels():
        tick.set_fontfamily('serif')
    for tick in barcode_plot.get_yticklabels():
        tick.set_fontfamily('serif')
    #barcode_plot.axes.xaxis.set_visible(False)
    #barcode_plot.axes.yaxis.set_visible(False)
    # Construct legend manually because GUDHI default places it awkwardly.
    dimensions = list(set(item[0] for item in clean_intervals))
    colormap = plot.cm.Set1.colors
    barcode_plot.legend(
       handles=[
           mpatches.Patch(color=colormap[dim],
           label=str(dim)) for dim in dimensions],
       loc="lower left",fontsize='20')
    # frameon=False
    # barcode_plot.spines['top'].set_visible(False)
    # barcode_plot.spines['right'].set_visible(False)
    # barcode_plot.spines['bottom'].set_visible(False)
    # barcode_plot.spines['left'].set_visible(False)
    if args.write:
        barcode_file = args.output_file + '_barcode.pdf'
        print('Saving barcode plot to ' + barcode_file)
        plot.savefig(barcode_file)

if args.diagram:
    print('...Generating persistence diagram.')
    clean_intervals = du.remove_inf_bars(intervals)
    diagram_plot = gudhi.plot_persistence_diagram(
        clean_intervals,
        legend=True,
        max_intervals=0)
    # Set formatting choices for diagrams.
    diagram_plot.set_title(title, fontsize='30', fontfamily='serif')
    diagram_plot.tick_params(labelsize='15')
    for tick in diagram_plot.get_xticklabels():
        tick.set_fontfamily('serif')
    for tick in diagram_plot.get_yticklabels():
        tick.set_fontfamily('serif')
    diagram_plot.set_xlabel("Birth", fontsize='16', fontfamily='serif')
    diagram_plot.set_ylabel("Death", fontsize='16', fontfamily='serif')
    # Manually set axis tick labels because we don't need the GUDHI
    # default three decimal places of precision.
    # See sourcecode for GUDHI.persistence_graphical_tools for details.
    (min_birth, max_death) = mb_md(clean_intervals,0.0)
    inf_delta = 0.1
    delta = (max_death - min_birth) * inf_delta
    infinity = max_death + 2*delta
    axis_end = max_death + delta / 2
    yt = diagram_plot.get_yticks()
    yt = yt[np.where(yt < axis_end)]
    yt = np.append(yt, infinity)
    ytl = ["%d" % e for e in yt]
    ytl[-1] = r'$\infty$'
    diagram_plot.set_yticklabels(ytl)
    if args.write:
        diagram_file = args.output_file + '_diagram.pdf'
        print('Saving diagram plot to ' + diagram_file)
        plot.savefig(diagram_file)

if args.bars or args.diagram:
    plot.show()

print('######################################################################')
