# lagrangian-betweenness

## Overview

The repository `lagrangian-betweenness` contains the codes designed to compute Lagrangian betweenness. When using this code, please acknowledge the authors by citing  [Ser-Giacomi et al. (2021)](#references).





## Instructions

1 - Change line 51 and 52 in `betw_computation.cpp` with the actual path to the folder named "Outputs", which is found inside code_Betweenness folder

2 - Do the same thing in the file `parameters.cpp` at line 322

3 - Compile the code with a C++ compiler (this can be done with the Makefile that will produce the executable `betw_computation.out`) and run the resulting executable. 



## Notes
The `betw_computation.out` function needs 2 files to be run:

- `parameters00000001.dat` is the file in which you insert the parameters: (1) `vel_product`: which hydrodynamical field you want to use (if you have more than one); 
(2) `verbose`: 1 if you want the code to be verbose, 0 otherwise; (3) `start_date`: the starting date of the advection; (4) `tau`: the advection time period; (5) `delta0_ftle`: the initial separation used to compute the ftle, which are used to compute the betweenness field; (6) `rand_idx`: this parameter needs to be used only if you run the code more than once at the same time. `rand_idx` determines the suffix of the `parameters00000001.dat`. In this case, `rand_idx=1`. If `rand_idx=2`, then the parameter file should be named `parameters00000002.dat`

- `Glob_Ekman_dt_00000001.trac` is the file that contain the list of points for which we want to compute the betweenness. It is formed by two columns, the first containing the latitude of the points, the second the longitude.
The name of this file must (i) begin with the name of the velocity field product used, which in this case is specified at line 319 of `parameters.cpp`; (ii) end with `rand_idx` used.


## References

[[Ser-Giacomi et al. 2021]](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.103.042309) Ser-Giacomi, E., et al. Lagrangian betweenness: a measure of bottlenecks in dynamical systems with oceanographic examples. *Nature Communications* (2021)



