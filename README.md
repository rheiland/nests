# nests

This repository attempts to create a model of (potentially) spreading "nests" of cells. There are two cell types: dividing and differentiated. 

The code in this repo was created using a branch of PhysiCell that provided cell differentiation using probabilities:
https://github.com/MathCancer/PhysiCell/tree/differentiation-release

The primary files of interest to the modeler are `/custom_modules/custom_nests.{h,cpp}`, which is where one:
1) creates the initial cell population - their cell types and locations,
2) creates the initial microenvironment (whatever substrates may exist),
3) creates most of the model-specific rules and behaviors.
