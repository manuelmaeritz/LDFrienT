# LDFrienT
LDFrienT was created within the scope of my master's project at the institute of applied physics at the University of Tübingen in 2020.
This project, on the one hand, provides a general framework for lattice density functional theory (LDFT). This framework comes in the form of an abstract base class, the class `LdftModel` in the module **ldft_model.py**. It provides some general functions to LDFT (such as Picard iterations, calculation of surface properties, etc.) and, more importantly, a nice structure for your LDFT-project. To take advatage this framework, simply create a class for your LDFT-model, which inherits from `LdftModel` and implement the abstract methodes. Those are specific to your model (e.g. the calculation of the free energy). The purpose of this framework was to make the code clean and uniform for every project of this type.

On the other hand, LDFrienT comes with four readily implemented LDFT-models which make use of the `LdftModel`-framework. Those implemented models are the
- two-dimensional lattice gas in the mean field approximation (the class `LG2dMf` in the module **lg_2d_mf.py**)
- two-dimensional lattice gas in the highlander approximation (the class `LG2dAOHighl` in the module **lg_2d_highl.py**)
- three-dimensional lattice gas in the mean field approximation (the class `LG3dMf` in the module **lg_3d_mf.py**)
- three-dimensional lattice gas in the highlander approximation (the class `LG3dAOHighl` in the module **lg_3d_highl.py**)

Moreover, the **generator.py**-skript provides some usefull functions to generate series of density profiles (samples) of a certain LDFT-model. Some of these functions are customized to the above classes (`LG2dMf`, `LG2dAOHighl`, `LG3dMf`, `LG3dAOHighl`), though may easily be extended to your LDFT-model.

## Requirements and dependencies
This project requires phython 3 and makes use of the following packages
- numpy
- matplotlib
- os
- picle
- functools
- abc
- scipy.optimize
- sys

If one wants to also take advantage of the notebooks, Jupyter notebook is required on your machine.

## Exaples and Tests
An easy example of how to use the `LdftModel`-framework properly is provided by the two classes `LG2dMf` and `LG3dMf`. A more complex example is given by `LG2dAOHighl` and `LG3dAOHighl`.

The notebooks **LG2dAOHighl-testsheet.ipynb** and **LG2dMeanField-testsheed.ipynb** were originaly intended for testing the classes `LG2dAOHighl` and `LG2dMf`. They might, however, also be used as a little tutorial to those classes.

The notebook **generator-testsheet.ipynb** was for testing the functions provided in **generator.py**, though serve as a tutorial to this module as well. Moreover, there is a notebook **generator-interface.ipynb**, which can be used, if one does not want to call the functions in **generator.py** from the terminal.

## Content and file structure
- **code**.............................................................................This folder holds all the code
  - **ldft_classes_v2** 
    - ldft_model.py................................................. contains ``LdftModel``-class (abstract class serving as LDFT framework)
    - lg_2d_highl.py................................................ contains ``LG2dAOHighl``-class (two-dimensional lattice gas in Highlander approximation)
    - lg_2d_mf.py................................................... contains ``LG2dMf``-class (two-dimensional lattice gas in mean field approximation)
    - lg_3d_higl.py................................................. contains ``LG3dAOHighl``-class (three-dimensional lattice gas in Highlander approximation)
    - lg_3d_mf.py................................................... contains ``LG3dMf``-class (three-dimensional lattice gas in mean field approximation)
  - LG2dAOHighl-testsheet.ipynb.............................. Test or tutorial for the ``LG2dAOHighl``-class
  - LG2dMeanField-testsheet.ipynb........................... Test or tutorial for the ``LG2dMeanField``-class
  - generator-interface.ipynb...................................... "Interface"-notebook for ``generator.py``-modul
  - generator-testsheet.ipynb..................................... Test or tutorial for the ``generator.py``-modul
  - generator.py.......................................................... contains functions to produce samples
- **samples**...................................................................... some exemplary samples produce by the code
- .gitignore
- License
- README_GitHub......................................................... This read-me file

## Documentation

