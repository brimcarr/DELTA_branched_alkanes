.. DeltaPersistence documentation master file, created by
   sphinx-quickstart on Mon Jun 22 19:17:45 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Delta Persistence Documentation
============================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   alkaneenergy
   sublevelpersistence  
   datautils

This package implements the persistent homology code from the NSF-DELTA project.

Overview
========

The project "Descriptors of Energy Landscapes using Topological Data Analysis" (DELTA) seeks to use topological data analysis to study the energy landscapes of chemical systems.
This package contains sourcecode used by the team at Colorado State University for analysis of energy landscapes of alkane molecules.
The base persistent homology computations are performed by the GUDHI package.
This package currently contains three modules which supplement the base GUDHI modules.

* Alkane Energy - provides methods for computing the energy landscapes of n-alkane series.
* Sublevel Persistence - provides general methods for computing sublevel set persistent homology.
* Data Utilities - provides methods for converting data to formats usable by GUDHI and for working with GUDHI output.

Getting started
===============

For users of delta-topology.org, see the section :ref:`deltagateway`.
For Delta project members wanting to work from the terminal on Stampede, see :ref:`stampede`.
If you want to work on your local machine, follow the below instructions.

First, the GUDHI library and several other scientific packages for python need to be installed.
This is easiest if using conda/anaconda.
Here is a minimal installation procedure for GUDHI:

Create a clean python 3.8 environment, here called "deltagudhi", though the name is flexible::

    conda create --name deltagudhi python=3.8

Make that the active environment::

    conda activate deltagudhi

List all the installed packages::

    conda list

At this point there should be basically nothing there. Next, install numpy and matplotlib, which are strict required for both GUDHI and this package::

    conda install numpy
    conda install matplotlib

Answer yes to any prompts to install other dependencies.
Install GUDHI from the conda-forge repository::

    conda install -c conda-forge gudhi

To check that everything works, start python::

    python3

Once in python, try importing gudhi::

    import gudhi

If everything went correctly the cursor should advance to the next line in the python interpreter with no output. To double-check, run::

    print(gudhi.__debug_info__)

The output should look something like this::

    Pybind11 version 2.5.0
    Python version 3.8.2
    Cython version 0.29.18
    Numpy version 1.14.6
    Eigen3 version 3.3.7
    CGAL header only version 5.00.1.100
    GMP_LIBRARIES = /opt/anaconda2/envs/marktest/lib/libgmp.dylib
    GMPXX_LIBRARIES = /opt/anaconda2/envs/marktest/lib/libgmpxx.dylib
    MPFR_LIBRARIES = /opt/anaconda2/envs/marktest/lib/libmpfr.dylib

At this point GUDHI is installed.
To use this package, download/clone the repository and navigate to the folder `deltapersistence/` (not the subfolder with the same name).
Run some basic tests to verify that everything works.
In detail, at the command prompt, run::

    python3 persistence_script.py --carbons 5 --disp

This should produce the persistence intervals for the energy landscape of pentane (the n-alkane with 5 carbon atoms).

Run::

    python3 persistence_script.py -f tests/cube.txt -t cube --bars

This should generate the bars for pentane, now using the detailed formula which accounts for non-bonded interactions. (This is the energy landscape saved as a cubical complex in `cube.txt`.)

Run::

    python3 persistence_script.py -f tests/mesh.npy -t mesh -pd 111 --diagram

This should generate the persistence diagram for hexane from the energy landscape saved as `mesh.npy`. It will take a few moments.

Run::

    python3 persistence_script.py -f tests/coords.dat -e tests/energy.dat -t sample --bars

to generate the barcodes for the sublevelsets of the function :math:`f(x,y) = - (x-.5)^2 - (y-.5)^2` on the square :math:`I = [0,1]\times[0,1]`.
The file `coords.dat` contains a random sample of points from *I* and `energy.dat` contains the value of *f* on those samples. 
Sublevelset persistence is computed by forming the Delaunay triangulation of the sampled coordinates, filtered by the lower-star filtration for *f*.

Finally, run::

    python3 persistence_script.py --help

This shows all available options.
Experiment with some of them!

Assuming that all goes well, the package can be installed for system-wide access.
To do this, from within the folder `deltapersistence` run::

    pip install .

Doing so allows you to import `deltapersistence` and its provided modules while working in any folder on your system.
It also installs the script `persistence_script.py` as a command-line utility which can be run from anywhere by calling `persistence_script`.

If you need capabilities beyond those accessible from `persistence_script.py` you can look at the detailed documentation of the modules :doc:`alkaneenergy`, :doc:`sublevelpersistence`, and :doc:`datautils`.

.. _deltagateway:
Delta Gateway
-------------

At https://delta-topology.org/, log in. Now at the DELTA Gateway dashboard, click on the project "DeltaPersistence". 

**Experiment 1: Pentane.**

This experiment computes the sublevelset persistent homology of the n-alkanes.
You get to choose the number of carbons in the chain as input --- 4 for butane, 5 for pentane, etc!
Under execution choice, select "n-alkane input".
Under alkane chain length, type in "5" for pentane.
Click the boxes to "Show Barcode", "Show Diagram", and "List Intervals" (these should automatically be selected).
Click "Save and Launch".
Once the experiment ends, under "Delta Persistence Output" click on "intervals_barcode.pdf".
These are the sublevelset persistent homology intervals, which should match those at 
https://gitlab.com/thrust-2/thrust2/-/blob/master/notes/V_pentane_periodic_barcodes.png,
which are explained in Figures 5 and 6 of the Alkane Notes:
https://gitlab.com/thrust-2/thrust2/-/blob/master/notes/AlkanePaperDraft_compiled.pdf.
Also under "Delta Persistence Output", click on "intervals_diagram.pdf" to see the persistent homology instead displayed as a persistence diagram.
The file "intervals.pers" is text output which you can also see by clicking on "DeltaPersistence.stdout" a few rows down.

*Advanced comments:* You can try changing the number of carbons in the alkane length, starting from 4 carbons (butane).
As the number of carbons increases, to reduce computation time you may want to reduce the resolution of the grid we place in the energy landscape; this is an optional parameter which one can add in the "advanced arguments".

**Experiment 2: Mesh data.**
This experiment computes the sublevelset persistent homology of real-valued function defined on a mesh, i.e. a 3-dimensional array / matrix.
The particular mesh that we use for this example stores the energy landscape for hexane, but you can instead input your own data.
Under execution choice, select "Mesh data input".
Download the file https://gitlab.com/thrust-2/thrust2/-/blob/master/software/deltapersistence/tests/mesh.npy and add it using the "Drop files here or browse" button on the DELTA Gateway.
Click the boxes to "Show Barcode", "Show Diagram", and "List Intervals" (these should automatically be selected).
Under "Advanced Arguments", type in the flag "-pd 111".
This tells GUDHI to add periodic boundary conditions in all three dimensions.
If you instead only wanted periodic boundary conditions in the first and third dimensions, you would type "-pd 101".
Click 'Save and Launch"
Once the experiment ends, under "Delta Persistence Output" click on "intervals_barcode.pdf".
These are the sublevelset persistent homology intervals, which should match those at 
https://gitlab.com/thrust-2/thrust2/-/blob/master/notes/V_hexane_periodic_barcodes.png,
which are explained in Figure 6 of the Alkane Notes:
https://gitlab.com/thrust-2/thrust2/-/blob/master/notes/AlkanePaperDraft_compiled.pdf.
Also under "Delta Persistence Output", click on "intervals_diagram.pdf" to see the persistent homology instead displayed as a persistence diagram.
The file "intervals.pers" is text output which you can also see by clicking on "DeltaPersistence.stdout" a few rows down.

**Experiment 3: Sample data input.**
This experiment computes the sublevelset persistent homology of real-valued function defined on arbitrary sample of points in some domain,
perhaps produced by molecular dynamics simulations.
The first input will be the coordinates of the data points.
The second coordinate will be the energy associated to each data point.
The script computes the Delaunay complex of the data points, as shown in Figure 5(right) of
https://gitlab.com/thrust-2/thrust2/-/blob/master/notes/Thrust2Notes_compiled.pdf
It then computes the sublevelset persistent homology based on a linear interpolation of the energy values.
Under execution choice, select "Sample data input".
Download the file https://gitlab.com/thrust-2/thrust2/-/blob/master/software/deltapersistence/tests/coords.dat and add it using the "Drop files here or browse" button on the DELTA Gateway under "Sample Data Coordinates".
Download the file https://gitlab.com/thrust-2/thrust2/-/blob/master/software/deltapersistence/tests/energy.dat and add it using the "Drop files here or browse" button on the DELTA Gateway under "Sample Data Energy".
Click the boxes to "Show Barcode", "Show Diagram", and "List Intervals" (these should automatically be selected).
Click 'Save and Launch"
Once the experiment ends, under "Delta Persistence Output" click on "intervals_barcode.pdf".
These are the sublevelset persistent homology intervals.
Also under "Delta Persistence Output", click on "intervals_diagram.pdf" to see the persistent homology instead displayed as a persistence diagram.
The file "intervals.pers" is text output which you can also see by clicking on "DeltaPersistence.stdout" a few rows down.



.. _stampede:
Stampede
--------

This package is installed on Stampede.
To access it, use the Delta Topology login to stampede2 and activate the conda environment `deltagudhi`.
All modules can then be imported into any python project, or `persistence_script.py` can be run.
