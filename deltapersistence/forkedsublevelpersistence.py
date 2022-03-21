"""Tools for computing sublevelset persistence

The GUDHI library provides methods for computing persistence of general
filtered complexes (cubical and simplicial). This module provides tools
for determining the sublevelset filtration of functions (equivalently,
the lower-star filtration of a complex) in compatible formats.


    Copyright (C) 2020 Joshua Mirth, modified by Brittany Story 2022

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

* _array_neighbors
* cutoff
* lower_star_filter
* mesh_sublevel_persistence
* sublevel_complex

"""

# General libraries.
import numpy as np

# Third-party libraries.
import gudhi

def mesh_sublevel_persistence(fnx,periodic_dims=[],field=2):
    """Compute the sublevelset persistence of a function defined on a
    mesh.

    Constructs a filtered cubical complex from a function defined on a
    cubical mesh, then computes the sublevelset persistence intervals
    for that function. Cubical here can be a mesh defined on any
    rectangular grid in Euclidean space, possibly with periodic boundary
    conditions.

    Parameters
    ----------
    fnx : numpy ndarray
        Array of function values, assumed to live on an evenly-spaced
        mesh of a Euclidean rectangle.
    periodic_dims : list of bools, optional
        Dimensions on which to assume periodic boundary conditions.
        Default is empty, that is, no periodic boundaries. Length of the
        list must be equal to the dimension of fnx.
    field : int, optional
        Field in which to compute homology. Must be a prime number.
        Default is to use Z_2.

    Returns
    -------
    intervals : list
        List of persistence intervals. Each interval is a tuple of the
        form `(dim, (birth, death))` where `dim` is the integer
        corresponding to the homological dimension and `birth` and
        `death` are floats giving the birth and death times.

    """

    cplx = sublevel_complex(fnx)
    if fnx.ndim == len(periodic_dims):
        cubical_cplx = gudhi.PeriodicCubicalComplex(
            top_dimensional_cells=cplx,
            periodic_dimensions=periodic_dims)
    else:
        cubical_cplx = gudhi.CubicalComplex(top_dimensional_cells=cplx)
        print('Assuming that no dimensions are periodic.')
    intervals = cubical_cplx.persistence(
        min_persistence=0,
        homology_coeff_field=field)
    return intervals

def lower_star_delaunay(coords,fnx,periodic=False):
    """Constructs lower-star filtered Delaunay complex from a function.

    Parameters
    ----------
    coords : numpy array
        List of coordinates in n-dimensional Euclidean space. Each row
        corresponds to one coordinate.
    fnx : numpy array
        List of function values at the coordinates.

    Returns
    -------
    lower_star_complex : GUDHI SimplexTree
        Delaunay complex of the coordinates with lower-star filtration.

    Raises
    ------
    ValueError: if number function array and coordinate array do not
        have the same number of rows.

    """

    if len(coords) != len(fnx):
        raise ValueError('Number of data points must be equal to number of ' +
            'function values.')
    if coords.ndim == 1:
        lower_star_complex = one_dimensional_delaunay(coords,periodic)
    elif periodic:
        # TODO: allow setting different values of epsilon here.
        lower_star_complex = periodic_delaunay(coords,1.0)
    else:
        cplx = gudhi.AlphaComplex(points=coords)
        lower_star_complex = cplx.create_simplex_tree()
    lower_star_filter(lower_star_complex,fnx)
    return lower_star_complex

def periodic_delaunay(coords,epsilon,ub=2*np.pi,lb=0,verbose=False):
    """Constructs an (approximate) Delaunay complex on a set of points
    with periodic boundary conditions.

    The points are assumed to lie in a (hyper)cube. To compute the
    Delaunay triangulation, the coordinates are tiled and the Euclidean
    Delaunay triangulation of the tiled data is computed. The simplices
    of this triangulation which lie near the boundary of the cube are
    then projected down onto the original dataset.

    Parameters
    ----------
    coords : numpy array
        List of coordinates in d-dimensional Euclidean space. Each row
        corresponds to one coordinate. Dimension (d) must be 2 or 3.
    epsilon : float
        Points within epsilon of the original cube will be projected
        when accounting for periodic boundaries. Epsilon should be
        positive.
    ub : float
        Upper bound for cube. Default is 2pi.
    lb : float
        Lower bound for cube. Default is 0.
    verbose : bool
        When true, print information about computations. Default is
        false.

    Returns
    -------
    delaunay_complex : GUDHI SimplexTree
        Periodic Delaunay complex of the input data.

    """

    if coords.ndim != 2:
        raise ValueError('Input data must be a two-dimensional numpy array ' +
            ' of coordinate values.')
    n = coords.shape[0]     # number of data points
    dim = coords.shape[1]
    duplicate_coords = np.zeros(((3**dim)*n,dim))
    width = ub - lb
    if dim == 2:
        shift_pattern = width*np.array(
            [[0,0],[-1,-1],[-1,0],[-1,1],[0,-1],[0,1],[1,-1],[1,0],[1,1]]
        )
    elif dim == 3:
        shift_pattern = width*np.array([
            [0,0,0],[-1,-1,-1],[-1,-1,0],[-1,-1,1],[-1,0,-1],[-1,0,0],[-1,0,1],
            [-1,1,-1],[-1,1,0],[-1,1,1],[0,-1,-1],[0,-1,0],[0,-1,1],[0,0,-1],
            [0,0,1],[0,1,-1],[0,1,0],[0,1,1],[1,-1,-1],[1,-1,0],[1,-1,1],
            [1,0,-1],[1,0,0],[1,0,1],[1,1,-1],[1,1,0],[1,1,1]])
    else:
        raise NotImplementedError('Periodic Delaunay only implemented for ' +
            ' dimensions two and three.')
    j = 0
    for p in shift_pattern:
        duplicate_coords[n*j:n*(j+1),:] = coords+p
        j += 1
    cplx = gudhi.AlphaComplex(points=duplicate_coords)
    simplex_tree = cplx.create_simplex_tree(default_filtration_value=True)
    if verbose:
        print('Original complex has ' + str(simplex_tree.num_vertices())
            + ' vertices and ' + str(simplex_tree.num_simplices()) +
            ' simplices.')
    simplices = [sigma for sigma in simplex_tree.get_simplices()]
    for sigma in simplices:
        pts = [cplx.get_point(x) for x in sigma[0]]
        if np.max(pts) > ub+epsilon or np.min(pts) < lb-epsilon:
            if verbose:
                print('Removing simplex ' + str(sigma[0]))
            simplex_tree.remove_maximal_simplex(sigma[0])
        elif np.max(pts) > ub or np.min(pts) < lb:
            if verbose:
                print('Projecting simplex ' + str(sigma[0]) + ' to ' +
                    str(np.mod(sigma[0],n)) + '.')
            simplex_tree.remove_maximal_simplex(sigma[0])
            t = simplex_tree.insert(np.mod(sigma[0],n),0)
            if verbose:
                if t:
                    print('Inserted simplex' + str(np.mod(sigma[0],n)))
                else:
                    print('Simplex ' + str(np.mod(sigma[0],n)) +
                        ' already exists.')
    if verbose:
        print('Updated complex has ' + str(simplex_tree.num_vertices()) +
            ' vertices and '+str(simplex_tree.num_simplices())+' simplices.')
    for sigma in simplex_tree.get_simplices():
        simplex_tree.assign_filtration(sigma[0],0)
    pers = simplex_tree.persistence(persistence_dim_max=True)
    if verbose:
        print('Homology of resulting complex: ' + str(pers))
    if dim == 2:
        expected_pers = [(2, (0.0, np.inf)), (1, (0.0, np.inf)),
            (1, (0.0, np.inf)), (0,(0.0, np.inf))]
    elif dim == 3:
        expected_pers = [(3, (0.0, np.inf)), (2, (0.0, np.inf)),
            (2, (0.0, np.inf)), (2, (0.0, np.inf)), (1, (0.0, np.inf)),
            (1, (0.0, np.inf)), (1, (0.0, np.inf)), (0, (0.0, np.inf))]
    if pers != expected_pers:
        print('WARNING: homology of complex is not equal to homology of the ' +
            ' torus. Consider trying a different epsilon value.')
    return simplex_tree

def lower_star_filter(simplex_tree,fnx):
    """Give a simplicial complex the lower-star filtration.

    The lower-star filtration assigns to each simplex the maximum
    filtration value of its vertices. Given any simplicial complex as a
    `simplex_tree`, and a function on its vertices, this function
    updates the filtration values to be the lower-star.

    Parameters
    ----------
    simplex_tree : GUDHI SimplexTree
        A filtered simplicial complex as a GUDHI simplex tree.
    fnx : list of floats
        List of function values on the vertices of simplex_tree.

    Notes
    -----
    SimplexTree is a mutable data type. Applying this function updates
    the value of the supplied object. GUDHI doe not provide a `copy()`
    method for objects of type `SimplexTree` as of version 3.2.0.

    """

    if len(fnx) != simplex_tree.num_vertices():
        raise ValueError('Number of function values does not match ' +
            'number of vertices.')
    # Erase current filtration values by setting all to -inf.
    for v in simplex_tree.get_simplices():
        simplex_tree.assign_filtration(v[0],-np.inf)
    # Update with filtration values from function.
    for i in range(0,len(fnx)):
        simplex_tree.assign_filtration([i],fnx[i])
    simplex_tree.make_filtration_non_decreasing()

def one_dimensional_delaunay(coords,periodic=True):
    """Constructs the trivial triangulation of a set of points on the
    real line.

    Parameters
    ----------
    coords : numpy array
        One dimensional array or list of coordinates on the real line.
    periodic : boolean, optional
        If true, make use periodic boundary conditions.

    Returns
    -------
    cplx : GUDHI SimplexTree
        Delaunay complex of the input data with all filtration values
        set to -infinity.

    """

    vertices = np.argsort(coords)
    cplx = gudhi.SimplexTree()
    for v in vertices:
        tmp = cplx.insert([v],-np.inf)
    for i in range(0,len(vertices)-1):
        tmp = cplx.insert([vertices[i],vertices[i+1]],-np.inf)
    if periodic:
        tmp = cplx.insert([vertices[0],vertices[-1]],-np.inf)
    return cplx


def sublevel_complex(fnx, levelsize=0.1, perseus=False):
    """Construct a sublevel filtered cubical complex from a function.

    Converts from a mesh of function values to a mesh of filtration
    times at which a function value occurs. Can be passed directly to
    GUDHI with gudhi.cubical_complex(top_dimensional_cells=cmplx).

    Parameters
    ----------
    fnx : numpy ndarray
        Values of a function on a coordinate mesh.
    levelsize : float, optional
        Resolution at which to filter. Default is 0.1 which is
        reasonable for most uses. Changing this value does not
        typically change performance of persistent homology
        computations.
    perseus : bool, optional.
        If true convert all filtration values to integers. Default is
        False. GUDHI can compute persistent homology with float-valued
        filtrations, while Perseus requires integer filtration values.

    Returns
    -------
    cmplx : numpy ndarray
        The filtration values of the cubical complex.

    """

    if perseus:
        fnx = np.int_(np.floor(fnx/levelsize))
    dimn = fnx.ndim
    cmplx = np.zeros(tuple(np.array(fnx.shape)-1))
    nbrs = _array_neighbors(dimn)
    it = np.nditer(cmplx, flags=['multi_index'])
    while not it.finished:
        mi = np.array(it.multi_index)
        cubemax = fnx[it.multi_index]
        for i in nbrs:
            cubemax = max(cubemax,fnx[tuple(mi+i)])
        cmplx[it.multi_index] = cubemax
        it.iternext()
    return cmplx

def _array_neighbors(n):
    """Creates all binary tuples of length n.

    The sublevel complex needs to check the value on a cube against all
    of its neighbors, that is, all elements of the array at position
    i+1 in some component(s).

    Parameters
    ----------
    n : int
        Dimension of array.

    Returns
    -------
    nbrs : list
        All binary tuples of length n.

    Examples
    --------
    >>> sublevelpersistence._array_neighbors(2)
    [[0, 0], [0, 1], [1, 0], [1, 1]]

    """

    dimn = str(n)
    nbrs = []
    for i in range(0,2**n):
        nbrs.append(list(format(i, '0'+dimn+'b')))
        nbrs[i] = [int(j) for j in nbrs[i]]
    return nbrs

def cutoff(fnx, threshold):
    """Replaces all function values above threshold with infinity.

    GUDHI considers a filtration value of inf to mean that the
    corresponding cell is not in the domain of the function. For non-
    bonded interactions, cells on the domain can have very large energy
    values which should be excluded.

    Parameters
    ----------
    fnx : numpy ndarray
        Array of function values.
    threshold : float
        Maximum function value to allow.

    Returns
    -------
    cutfnx : numpy ndarray
        function with large values replaced by inf.

    Examples
    --------
    >>> sublevelpersistence.cutoff(np.array([1,2,3]),2)
    array([ 1.,  2., inf])

    """

    infarray = fnx*np.inf
    cutfnx = np.where(fnx > threshold, infarray, fnx)
    return cutfnx

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
