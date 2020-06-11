# Cost/loss analysis for CME forecast properties

## Introduction

This repository provides an implementation of the Cost/Loss analysis in Python, as described by [Owens et al. (2020)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2020SW002507).

## Installation
The ``Cost/Loss`` analysis is written in Python 3.7.3 and requires ``numpy``, ``matplotlib``, ``pandas`` and ``astropy``. For convenience, all required data are included. OMNI data are available from https://omniweb.gsfc81.nasa.gov/. The updated near-Earth CME list is available from http://www.srl.caltech.edu/ACE/ASC/DATA/level3/icmetable2.htm.

The simplest way to work with ``CostLoss`` in ``conda`` is to create its own environment. With the anaconda prompt, in the root directory of ``CostLoss``, this can be done as:
```
>>conda env create -f environment.yml
>>conda activate costloss
``` 
Then the examples can be accessed through 
```
>>jupyter lab CostLoss.ipynb
```
## Contact
Please contact either [Mathew Owens](https://github.com/mathewjowens) or [Luke Barnard](https://github.com/lukebarnard). 

## Citation
Please cite this software as Owens, M. J., Lockwood, M., & Barnard, L. A. ( 2020). The value of CME arrival‐time forecasts for space weather mitigation. Space Weather, 18, e2020SW002507. https://doi.org/10.1029/2020SW002507
