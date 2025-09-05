# COBRA-k TODOs

## Short-term (0.0.5)

* X Compartmentalize kinetic value collection, making them work with COBRA-k models alone or also annotated SBMLs
* X Add references to annotated SBMLs
* Add "old" bottleneck function
* X Add in/out function in utilities, also as plot
* X Add IPOPT executable option
* Create COBRA-k template repository and reference in documentation/README

## Mid-term

* Add gallery in documentation for plot functions
* Add "CompactModel" class
* Consistent argument names
* Consistent argument orders
* Add rank function
* Add linear-fractional programming
* Community models
* Module-based MCS (StrainDesign integration?)
* MILP gap filling
* More and better error messages
* Simple MILP EFM function for very small networks
* Switch to StringZilla
* High code coverage
* prompt-toolkit interactive model view
* ASCII plots

## Long-term

* Performant EFM utility
* Look up usage of populate
* numba usage?
* REPL integration (in Spyder?)
* Formation energies as thermodynamic constraint alternative in Metabolite
* Cofactor swapping routines
* No-GIL support for Python ≥ 3.13
* Parallelize equilibrator handling and k_cat/k_m collection
