# gflex_load_SE

This code has been provided "as is", for reproducibility purposes (Solid Earth, "Glacial-isostatic-adjustment strain rate–stress paradox in the Western Alps and impact on active faults and seismicity"). 
We do not provide guarantees about the model's validity for any other use.


This code uses the gFlex code (Wickert, 2016) which solves the flexure equations in 3D, here from an analytical solution. The elastic thickness is therefore constant.
Was implemented:

- Reading a grid representing a spatially variable load
- Viscous lithospheric relaxation (see Turcotte and Schubert)
- Calculations of strain and stress tensors at any point


The modeling parameters can be modified in launcher-gflex.py. Except in specific cases, it is not necessary to change anything in gflex_load.py.
