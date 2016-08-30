# QuaHoG
### QCD library for Hadronic and Gauge measurements

Parallel (MPI + OpenMP) functions for hadronic matrix element calculations:
* Smearing of gauge links
* Preparing sources, including smeared sources and sequential sources
* Contracting propagators into correlation functions
* Reading and writing correlation functions into structured binary files, currently HDF5
* Does not provide solvers. Rather, relies on external solver libraries ([DDalphaAMG](https://github.com/sbacchio/DDalphaAMG) or [tmLQCD](https://github.com/etmc/tmLQCD)).

Main programs for driving QuaHoG functions can be found in the [hadstruct-progs](https://github.com/g-koutsou/hadstruct-progs) repo.
