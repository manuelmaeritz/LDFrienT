# The code folder
This folder contains all the code necessary for data production and some test notebooks.

## The ldft_classes_v2 folder
The realisation of a lattice model - lets say for example a two dimensional lattice gas in mean field approximation of size 40x40, temperature epsi=2.0 and average density 0.3 - is implemented as the instance of a class, here.
The model and the functional (in the above example 'lattice gas' and 'mean field') is determined by the class. The physical state of the system (size, temperature, density/chemical potential) is determined by the parameters of the specific instance.

### The `LdftModel`-class
The `LdftModel`-class is implemented in the module **ldft_model.py**. It is an abstract class, which means it can not be instantiated itself. I tried to keep it as general as possible and therefore might also be applied to other models beside the lattice gas.
It provides a data structure for lattice systems by requiring some instance-variables. Moreover some general functions, which do not depend on the specific model are implemented here (e.g. the picard iteration, calculation of surface properties).
For a complete list of instance-variables and provided methods please try the help function (`help(LdftModel)`). Note, however, that some of the listed methods are abstract, which means they are not yet implemented but need to be implemented by the inherited classes. To see which method is abstract or not you need to have a look in the source code.

### The `LG2dMf`-class
The `LG2dMf`-class is implemented in the module **lg_2d_mf.py**. It inherits from `LdftModel` and is the class for two-dimensional lattice gases in mean field approximation. Pleas use the `help(LG2dMf)` function to get an overview of the class. For the help regarding a specific function use `help(LG2dMf.function)`.

### The `LG3dMf`-class
The `LG3dMf`-class is implemented in the module **lg_3d_mf.py**. It inherits from `LdftModel` and is the class for three-dimensional lattice gases in mean field approximation. Pleas use the `help(LG3dMf)` function to get an overview of the class. For the help regarding a specific function use `help(LG3dMf.function)`.

### The `LG2dAOHighl`-class
The `LG2dAOHighl`-class is implemented in the module **lg_2d_highl.py**. It inherits from `LdftModel` and is the class for two-dimensional lattice gases in the Highlander-approximation. Pleas use the `help(LG2dAOHighl)` function to get an overview of the class. For the help regarding a specific function use `help(LG2dAOHighl.function)`.

### The `LG3dAOHighl`-class
The `LG3dAOHighl`-class is implemented in the module **lg_3d_highl.py**. It inherits from `LdftModel` and is the class for two-dimensional lattice gases in the Highlander-approximation. Pleas use the `help(LG3dAOHighl)` function to get an overview of the class. For the help regarding a specific function use `help(LG3dAOHighl.function)`.

## The generator module
The module **generator.py** provides some methods for generating instances of the classes above. This module can be either uses from the command line or you can load the functions of the module in your personal script and use them alone. For help how to use the script form command line, read the documentation in the header of the module. For help to a specific function call the help (`help(function)`).

## The generator-interface
The notebook **generator-interface.ipynb** provides sort of an interface for the functions defined in the **generator.py** module.

## The generator-testsheet
The notebook **generator-testsheet.ipynb** was originally intended to test the functions of the **generator.py** module. It, however, may also be used as a little tutorial for the **generator.py** script.

## The LG2dAOHighl-testsheet
The notebook **LG2dAOHighl-testsheet.ipynb** was originally used for testing the methods provided by the `LG2dAOHighl`-class. It, however, may also be used as a little tutorial for the this class.

## The LG2dMeanField-testsheet
The notebook **LG2dMeanField-testsheet.ipynb** was originally used for testing the methods provided by the `LG2dMeanField`-class. It, however, may also be used as a little tutorial for the this class.
