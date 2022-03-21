"""Utilities for handling GUDHI and DELTA data types.

Converts internal data types used by GUDHI and provides methods for
reading and writing different data formats used by the DELTA team. In
particular, the following standard formats are used:

=========================== ================    ================    =============
Semantic Type               Python data type    File Format         Specification
=========================== ================    ================    =============
Energy landscape (sampled)  Numpy array         Plaintext (`.dat`)  (Tab-separated values)
Energy landscape (meshed)   Numpy ndarray       Numpy (`.npy`)      https://numpy.org/doc/stable/reference/generated/numpy.lib.format.html#module-numpy.lib.format
Sublevelset cubical complex Numpy ndarray       Perseus (`.txt`)    https://gudhi.inria.fr/python/latest/fileformats.html#perseus
Persistence intervals       List of tuples      GUDHI (`.pers`)     https://gudhi.inria.fr/python/latest/fileformats.html#persistence-diagram
=========================== ================    ================    =============

    Copyright (C) 2020 Joshua Mirth and Johnathan Bush, modified by Brittany Story 2022

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

* array_to_barcodes
* downsample
* intervals_in_dimension
* is_perseus_file
* plot_persistence_barcode
* plot_persistence_diagram
* remove_inf_bars
* sequential_maxmin
* write_sublevels
* write_intervals

"""

# Standard packages
import numpy as np
import random
from pathlib import Path
from os import path

def write_intervals(intervals, filename, message='Persistence Intervals'):
    """Write the output persistence intervals to file.

    Takes a list of persistence intervals and writes them to a file in
    the persistence format specified by GUDHI. These can then be read
    and plotted using GUDHI's built-in methods.

    Parameters
    ----------
    filename : string
        Name of file to write intervals to. Should end with ".pers",
        though this is not enforced.
    intervals : list of tuples
        List of persistence intervals in the form (dim,(birth,death)).
    message : string, optional
        Message to print as the first line of the file.

    """

    formatted_intervals = ''
    for i in intervals:
        formatted_intervals = formatted_intervals + str(i) + '\n'
    formatted_intervals = formatted_intervals.replace(')','')
    formatted_intervals = formatted_intervals.replace('(','')
    formatted_intervals = formatted_intervals.replace(',','')
    with open(filename, 'w') as filehandle:
        filehandle.writelines('# ' + message + '\n')
        filehandle.writelines(item for item in formatted_intervals)

def array_to_barcodes(bar_array,dimn):
    """Convert an ndarray to a list of birth death times.

    When computing only one dimension of homology, GUDHI returns an
    array of birth/death times rather than a list of birth/death
    tuples. This method converts such arrays into valid birth/death
    lists so that standard interval methods can be used.

    Parameters
    ----------
    bar_array : n*2 numpy ndarray
        Array of birth death times.
    dimn : int
        Dimension for which these are barcodes.

    Returns
    -------
    intervals : list of tuples
        List of persistent homology intervals.

    Raises
    ------
    ValueError : if the input array is not size n*2 (required to be
        interpreted as barcodes).


    Notes
    -----
    This method should not be necessary. The GUDHI documentation
    specifies that an N*2 array of birth/death times is valid input to
    any of the persistence graphical tools (see
    http://gudhi.gforge.inria.fr/python/latest/persistence_graphical_tools_ref.html).
    However, these methods to not seem to actually work unless the input
    consists of a full persistence values list. Thus this method exists
    to essentially duplicate the private `_array_handler` method in the
    GUDHI graphical tools module.

    """

    if bar_array.shape[1] != 2:
        raise ValueError('Provided array not formatted as birth/death times.')
    intervals = []
    for r in  bar_array:
        intervals.append((dimn,tuple(r)))
    return intervals

def remove_inf_bars(persistence):
    """Remove any bars in the persistence which are born at infinity.

    When using the periodic cubical complex method it is possible for
    persistence intervals to be born at infinity. These give an error if
    passed to barcode or diagram plotting methods and need to be
    removed.

    Parameters
    ----------
    persistence : list of tuples
        A list of tuples of the form (dimn, (birth, death)).

    Returns
    -------
    intervals : list of persistence
        The input list, but with any entries where birth=infinity
        removed.

    """
    # Assumes the persistence is sorted by dimension from largest to
    # smallest (as it should be from all standard methods).
    intervals = persistence.copy()  # Do not treat the input as mutable.
    max_dim = intervals[0][0]
    for i in range (0,max_dim+1):
        try:
            intervals.remove((i,(np.inf,np.inf)))
            # TODO: probably better to pass this message as a return
            # value.
            print('Removed interval born at infinity in dimension ' +
                str(i))
        except ValueError:
            pass
    return intervals

def intervals_in_dimension(persistence,dimn):
    """Returns the persistence intervals in the given dimension.

    When GUDHI computes persistence it always computes in all
    dimensions, however sometimes it is desirable to only return the
    intervals in a single dimension (e.g. if there are a very large
    number of total intervals).

    Parameters
    ----------
    persistence : list of tuples
        A list of tuples of the form (dimn, (birth, death)).
    dimn : int
        Dimension of homology to return.

    Returns
    -------
    intervals : list of tuples
        The list of intervals in the form (dimn, (birth, death)) where
        dimn is restricted to the specified dimension.

    Notes
    -----
    This function is not strictly necessary since most GUDHI objects
    possessing persistence have a `persistence_intevals_in_dimension()`
    method. However, sometimes it is convenient to work with the
    already-returned list of intervals rather than the complex object.

    """

    intervals = []
    for i in persistence:
        if i[0] == dimn:
            intervals.append(i)
    return intervals

def downsample(coords_filepath, energy_filepath, n):
    """Downsamples coordinate and energy data.

    Downsamples (uniformly) two data files, one containing function
    coordinates and one containing function (energy) values, and writes
    the result to a new pair of files, maintaining the ordering and
    pairing of coordinates with energy values.

    Parameters
    ----------
    coords_filepath : string
        Filepath to coordinate data.
    energy_filepath : string
        Filepath to energy data.
    n : integer
        Number of entries to be randomly sampled from both lists.

    Notes
    -----
    Output is a pair of files in the current working directory, with
    file names inherited from the input files, prepending 'downsampled'
    and the number of resulting data points.

    """

    coords_count = 0
    energy_count = 0
    with open(coords_filepath, 'r') as f:
        for line in f:
            coords_count += 1
    with open(energy_filepath, 'r') as f:
        for line in f:
            energy_count += 1
    if coords_count != energy_count:
        raise ValueError('Input files must be the same size.')
    random_lines = sorted(random.sample(range(0, coords_count), n))
    print(random_lines)
    filecoords = open(coords_filepath, 'r')
    coordslines = filecoords.read().splitlines()
    fileenergy = open(energy_filepath, 'r')
    energylines = fileenergy.read().splitlines()
    coords_outfile = 'downsample_' + str(n) + '_' + Path(coords_filepath).name
    energy_outfile = 'downsample_' + str(n) + '_' + Path(energy_filepath).name
    with open(coords_outfile, 'w') as output_file:
        output_file.writelines(
            coordslines[line] + '\n' for line in random_lines
        )
    with open(energy_outfile, 'w') as output_file:
        output_file.writelines(
            energylines[line] + '\n' for line in random_lines
        )
    filecoords.close()
    fileenergy.close()

def sequential_maxmin(coords,subsample_size,seed=0):
    """Use a sequential maxmin algorithm to downsample input data.

    Parameters
    ----------
    coords : numpy float array
        Array of coordinates in d-dimensional space.
    subsample_size : integer
        Size of desired downsample.
    seed : integer, optional
        Index of initial point. Default is to use the first point in
        the data. Users may desire to choose this point randomly.

    Returns
    -------
    maxmin_indices : integer list
        List containing indices of the coordinates selected by maxmin.

    """

    maxmin_indices = [seed]
    distances = np.full(len(coords),np.inf)
    while len(maxmin_indices) < subsample_size:
        pt = coords[seed]
        distances = np.minimum(np.sum((coords - pt)**2,axis=1), distances)
        seed = distances.argmax()
        maxmin_indices.append(seed)
    return maxmin_indices

########################################################################
# WARNING: The remaining code is copied from the GUDHI project with
# small modifications. It should be removed and replaced with calls to
# the actual GUDHI library. This is currently necessary because GUDHI
# does not properly handle not having TeX installed, and because we want
# to change the formatting of barcodes in certain ways which are not
# easily done from outside the plot_barcodes function.
#
# Current code from GUDHI 3.2.0.
########################################################################

def __min_birth_max_death(persistence, band=0.0):
    """This function returns (min_birth, max_death) from the persistence.

    :param persistence: The persistence to plot.
    :type persistence: list of tuples(dimension, tuple(birth, death)).
    :param band: band
    :type band: float.
    :returns: (float, float) -- (min_birth, max_death).
    """
    # Look for minimum birth date and maximum death date for plot optimisation
    max_death = 0
    min_birth = persistence[0][1][0]
    for interval in reversed(persistence):
        if float(interval[1][1]) != float("inf"):
            if float(interval[1][1]) > max_death:
                max_death = float(interval[1][1])
        if float(interval[1][0]) > max_death:
            max_death = float(interval[1][0])
        if float(interval[1][0]) < min_birth:
            min_birth = float(interval[1][0])
    if band > 0.0:
        max_death += band
    return (min_birth, max_death)


def _array_handler(a):
    '''
    :param a: if array, assumes it is a (n x 2) np.array and return a
                persistence-compatible list (padding with 0), so that the
                plot can be performed seamlessly.
    '''
    if isinstance(a[0][1], np.float64) or isinstance(a[0][1], float):
        return [[0, x] for x in a]
    else:
        return a


def plot_persistence_barcode(
    persistence=[],
    persistence_file="",
    alpha=1.0,
    max_intervals=1000,
    max_barcodes=1000,
    inf_delta=0.1,
    legend=False,
    colormap=None,
    axes=None,
    fontsize=16,
):
    """This function plots the persistence bar code from persistence values list
    , a np.array of shape (N x 2) (representing a diagram
    in a single homology dimension),
    or from a `persistence diagram <fileformats.html#persistence-diagram>`_ file.

    :param persistence: Persistence intervals values list. Can be grouped by dimension or not.
    :type persistence: an array of (dimension, array of (birth, death)) or an array of (birth, death).
    :param persistence_file: A `persistence diagram <fileformats.html#persistence-diagram>`_ file style name
        (reset persistence if both are set).
    :type persistence_file: string
    :param alpha: barcode transparency value (0.0 transparent through 1.0
        opaque - default is 0.6).
    :type alpha: float.
    :param max_intervals: maximal number of intervals to display.
        Selected intervals are those with the longest life time. Set it
        to 0 to see all. Default value is 1000.
    :type max_intervals: int.
    :param inf_delta: Infinity is placed at :code:`((max_death - min_birth) x
        inf_delta)` above :code:`max_death` value. A reasonable value is
        between 0.05 and 0.5 - default is 0.1.
    :type inf_delta: float.
    :param legend: Display the dimension color legend (default is False).
    :type legend: boolean.
    :param colormap: A matplotlib-like qualitative colormaps. Default is None
        which means :code:`matplotlib.cm.Set1.colors`.
    :type colormap: tuple of colors (3-tuple of float between 0. and 1.).
    :param axes: A matplotlib-like subplot axes. If None, the plot is drawn on
        a new set of axes.
    :type axes: `matplotlib.axes.Axes`
    :param fontsize: Fontsize to use in axis.
    :type fontsize: int
    :returns: (`matplotlib.axes.Axes`): The axes on which the plot was drawn.
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib import rc
        plt.rc('text', usetex=False)
        plt.rc('font', family='serif')

        if persistence_file != "":
            if path.isfile(persistence_file):
                # Reset persistence
                persistence = []
                diag = read_persistence_intervals_grouped_by_dimension(
                    persistence_file=persistence_file
                )
                for key in diag.keys():
                    for persistence_interval in diag[key]:
                        persistence.append((key, persistence_interval))
            else:
                print("file " + persistence_file + " not found.")
                return None

        persistence = _array_handler(persistence)

        if max_barcodes != 1000:
            print("Deprecated parameter. It has been replaced by max_intervals")
            max_intervals = max_barcodes

        if max_intervals > 0 and max_intervals < len(persistence):
            # Sort by life time, then takes only the max_intervals elements
            persistence = sorted(
                persistence,
                key=lambda life_time: life_time[1][1] - life_time[1][0],
                reverse=True,
            )[:max_intervals]

        if colormap == None:
            colormap = plt.cm.Set1.colors
        if axes == None:
            fig, axes = plt.subplots(1, 1)

        persistence = sorted(persistence, key=lambda birth: birth[1][0])

        (min_birth, max_death) = __min_birth_max_death(persistence)
        ind = 0
        delta = (max_death - min_birth) * inf_delta
        # Replace infinity values with max_death + delta for bar code to be more
        # readable
        infinity = max_death + delta
        axis_start = min_birth - delta
        # Draw horizontal bars in loop
        for interval in reversed(persistence):
            if float(interval[1][1]) != float("inf"):
                # Finite death case
                axes.barh(
                    ind,
                    (interval[1][1] - interval[1][0]),
                    height=0.8,
                    left=interval[1][0],
                    alpha=alpha,
                    color=colormap[interval[0]],
                    linewidth=1,
                    align='edge',
                )
            else:
                # Infinite death case for diagram to be nicer
                axes.barh(
                    ind,
                    (infinity - interval[1][0]),
                    height=0.8,
                    left=interval[1][0],
                    alpha=alpha,
                    color=colormap[interval[0]],
                    linewidth=1,
                    align='edge',
                )
            ind = ind + 1

        if legend:
            dimensions = list(set(item[0] for item in persistence))
            axes.legend(
                handles=[
                    mpatches.Patch(color=colormap[dim], label=str(dim))
                    for dim in dimensions
                ],
                loc="lower right",
            )

        axes.set_title("Persistence barcode", fontsize=fontsize)

        # Ends plot on infinity value and starts a little bit before min_birth
        axes.axis([axis_start, infinity, 0, ind])
        return axes

    except ImportError:
        print("This function is not available, you may be missing matplotlib.")



def plot_persistence_diagram(
    persistence=[],
    persistence_file="",
    alpha=0.6,
    band=0.0,
    max_intervals=1000,
    max_plots=1000,
    inf_delta=0.1,
    legend=False,
    colormap=None,
    axes=None,
    fontsize=16,
    greyblock=True
):
    """This function plots the persistence diagram from persistence values
    list, a np.array of shape (N x 2) representing a diagram in a single
    homology dimension, or from a `persistence diagram <fileformats.html#persistence-diagram>`_ file`.

    :param persistence: Persistence intervals values list. Can be grouped by dimension or not.
    :type persistence: an array of (dimension, array of (birth, death)) or an array of (birth, death).
    :param persistence_file: A `persistence diagram <fileformats.html#persistence-diagram>`_ file style name
        (reset persistence if both are set).
    :type persistence_file: string
    :param alpha: plot transparency value (0.0 transparent through 1.0
        opaque - default is 0.6).
    :type alpha: float.
    :param band: band (not displayed if :math:`\leq` 0. - default is 0.)
    :type band: float.
    :param max_intervals: maximal number of intervals to display.
        Selected intervals are those with the longest life time. Set it
        to 0 to see all. Default value is 1000.
    :type max_intervals: int.
    :param inf_delta: Infinity is placed at :code:`((max_death - min_birth) x
        inf_delta)` above :code:`max_death` value. A reasonable value is
        between 0.05 and 0.5 - default is 0.1.
    :type inf_delta: float.
    :param legend: Display the dimension color legend (default is False).
    :type legend: boolean.
    :param colormap: A matplotlib-like qualitative colormaps. Default is None
        which means :code:`matplotlib.cm.Set1.colors`.
    :type colormap: tuple of colors (3-tuple of float between 0. and 1.).
    :param axes: A matplotlib-like subplot axes. If None, the plot is drawn on
        a new set of axes.
    :type axes: `matplotlib.axes.Axes`
    :param fontsize: Fontsize to use in axis.
    :type fontsize: int
    :param greyblock: if we want to plot a grey patch on the lower half plane for nicer rendering. Default True.
    :type greyblock: boolean
    :returns: (`matplotlib.axes.Axes`): The axes on which the plot was drawn.
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib import rc
        plt.rc('text', usetex=False)
        plt.rc('font', family='serif')

        if persistence_file != "":
            if path.isfile(persistence_file):
                # Reset persistence
                persistence = []
                diag = read_persistence_intervals_grouped_by_dimension(
                    persistence_file=persistence_file
                )
                for key in diag.keys():
                    for persistence_interval in diag[key]:
                        persistence.append((key, persistence_interval))
            else:
                print("file " + persistence_file + " not found.")
                return None

        persistence = _array_handler(persistence)

        if max_plots != 1000:
            print("Deprecated parameter. It has been replaced by max_intervals")
            max_intervals = max_plots

        if max_intervals > 0 and max_intervals < len(persistence):
            # Sort by life time, then takes only the max_intervals elements
            persistence = sorted(
                persistence,
                key=lambda life_time: life_time[1][1] - life_time[1][0],
                reverse=True,
            )[:max_intervals]

        if colormap == None:
            colormap = plt.cm.Set1.colors
        if axes == None:
            fig, axes = plt.subplots(1, 1)

        (min_birth, max_death) = __min_birth_max_death(persistence, band)
        delta = (max_death - min_birth) * inf_delta
        # Replace infinity values with max_death + delta for diagram to be more
        # readable
        infinity = max_death + 2*delta
        axis_end = max_death + delta / 2
        axis_start = min_birth - delta

        # bootstrap band
        if band > 0.0:
            x = np.linspace(axis_start, infinity, 1000)
            axes.fill_between(x, x, x + band, alpha=alpha, facecolor="red")
        # lower diag patch
        if greyblock:
            axes.add_patch(mpatches.Polygon([[axis_start, axis_start], [axis_end, axis_start], [axis_end, axis_end]], fill=True, color='lightgrey'))
        # Draw points in loop
        pts_at_infty = False  # Records presence of pts at infty
        for interval in reversed(persistence):
            if float(interval[1][1]) != float("inf"):
                # Finite death case
                axes.scatter(
                    interval[1][0],
                    interval[1][1],
                    alpha=alpha,
                    color=colormap[interval[0]],
                )
            else:
                pts_at_infty = True
                # Infinite death case for diagram to be nicer
                axes.scatter(
                    interval[1][0], infinity, alpha=alpha, color=colormap[interval[0]]
                )
        if pts_at_infty:
            # infinity line and text
            axes.plot([axis_start, axis_end], [axis_start, axis_end], linewidth=1.0, color="k")
            axes.plot([axis_start, axis_end], [infinity, infinity], linewidth=1.0, color="k", alpha=alpha)
            # Infinity label
            yt = axes.get_yticks()
            yt = yt[np.where(yt < axis_end)] # to avoid ploting ticklabel higher than infinity
            yt = np.append(yt, infinity)
            ytl = ["%d" % e for e in yt]
            ytl[-1] = r'$\infty$'
            axes.set_yticks(yt)
            axes.set_yticklabels(ytl)

        if legend:
            dimensions = list(set(item[0] for item in persistence))
            axes.legend(
                handles=[
                    mpatches.Patch(color=colormap[dim], label=str(dim))
                    for dim in dimensions
                ]
            )

        axes.set_xlabel("Birth", fontsize='20')
        axes.set_ylabel("Death", fontsize='20')
        axes.set_title("Persistence diagram", fontsize=fontsize)
        # Ends plot on infinity value and starts a little bit before min_birth
        axes.axis([axis_start, axis_end, axis_start, infinity + delta/2])
        return axes

    except ImportError:
        print("This function is not available, you may be missing matplotlib.")
