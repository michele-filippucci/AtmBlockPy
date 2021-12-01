# AtmBlockPy
Tools for atmospheric blocking data analysis written in Python.

Structure:

/lib :: contains classes:
	- BlockTools.py: toolkit for atmosphere blocking analysis
	- BlockPlots.py: file containing some functions for performing
      plots
          
    Function descriptions are found in the library script

main folder:: contains example of action performed with libraries.

This library is a set of tools for the analysis of atmospheric blocking in the northern hemisphere.
The index used for atm blocking diagnostic is described in "Davini et al. - 2012 - Bidimensional diagnostics, variability, and 
trends of northern hemisphere blocking". Some differences and features are added: the persistence and area criteria
are applied at the level of tracking. Tracking also allow the user to perform lagrangian analysis.

The requirements for this library are:
Python              3.8.10
xarray              0.18.2
numpy               1.20.3
scipy               1.6.3
tqdm                4.61.1
matplotlib          3.4.2
Cartopy             0.19.0.post1



