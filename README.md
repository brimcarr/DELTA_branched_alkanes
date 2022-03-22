# DELTA_branched_alkanes

Package implementing the persistent homology code from the NSF-DELTA project. This was modified to work with the added bond types.


## Outline

This python package contains six modules:
* Characterization for 32 - Gives the methods for distributing the critical points of an energy landscape for a molecule with internal bonds of type 32 into classes and calculates the number of critical points in each class. There is also a function that will characterize the sublevelset persistence barcode for energy landscapes for molecules consisting exclusively of internal bonds of type 3-2.
* Characterization for 32 22 - Same as above, but for molecules with internal bonds of types 3-2 and 2-2.
* Forked Alkane Energy - provides methods for computing the energy landscapes of the branched alkane series, using the formulae provided by Biswajit Sadhu.
* Forked Sublevel Persistence - provides general methods for computing sublevel set persistent homology (using GUDHI).
* Forked Data Utilities - provides a number of functions for converting data into formats GUDHI understands and manipulating GUDHI output for analysis.
* Equations of bonds of internal types 1-x-y-1: Code from Biswajit Sadhu to calculate the energy landscapes of all 6 molecules consisting exclusively of dihedral types 1-x-y-1.

There is also a command-line script, `persistence_script.py` which provides easy user access to a standard sublevelset persistence workflow.


## Contributing

We will use the "feature branch" workflow for git; see the useful links below.

Details of the development issues listed above will be added as issues in the Gitlab repository.

Periodically (when the code is "complete" up to some set of features) we will copy this to the main DELTA repository, but all development should happen here.


## Useful Links

* About documenting python code: https://realpython.com/documenting-python-code/
* Numpy-style docstrings: https://numpydoc.readthedocs.io/en/latest/format.html
* Python modules: https://www.tutorialspoint.com/python/python_modules.htm
* Structure of a python package: https://packaging.python.org/tutorials/packaging-projects/
* Tutorial on the git "feature branch" workflow: https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow
* Tutorial on general git collaboration: https://www.atlassian.com/git/tutorials/syncing
* Style guide for python: https://www.python.org/dev/peps/pep-0008/#function-and-variable-names

## Getting started

* Download and install GUDHI.
* Clone this repository.
* Navigate into the folder `deltapersistence/`.
* At the command prompt, run `python3 persistence_script.py --carbons 5 --disp`.
  This should produce the barcodes for pentane.
* Run `python3 persistence_script.py --bonds 1 1 0 0 0 0 --bars`.
  This should print the barcodes for an approximation of 2-methylpentane.
* Run `python3 persistence_script.py -f ../tests/cube.txt -t cube --bars`.
  This should generate the bars for pentane with non-bonded interactions which
  are saved as a cubical complex in the data file `cube.txt`.
* Run `python3 persistence_script.py -f ../tests/mesh.npy -t mesh --diagram`.
  This should generate the persistence diagram for hexane from the
  energy landscape saved as `mesh.npy`. It will take a few moments.
* Run `python3 persistence_script.py --help` to see all available
  options, then experiment with some of them!
  
## Documentation

For full documentation of the original code, see https://jrmirth.gitlab.io/deltapersistence/
For full documentation of the modified code, see https://github.com/brimcarr/DELTA_branched_alkanes/
This code accompanies the paper: Molecular configurations and persistence: Additive energies and branched alkanes by Henry Adams, Aurora Clark, Biswajit Sadhu, and Brittany Story.

