# About the samples file structure

The data is stored in the following file structure:
**model/size/temperature/**

**The model:** refers to the dft-model of the data (e.g. two dimensional Highlander model)

**The size:** refers to the system size (e.g. 64x64)

**The temperature** refers to the attraction strength times inverse temperature of the model (e.g. beta\*epsilon=2.0)

The filename usually starts with 'dens=' followed by a density. Optionally the density might be followed by a specification string in braces. The file-type is pickle. The saved files are instances of a LdftModel-class (see the README in the code-folder). 
