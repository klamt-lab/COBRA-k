# COBRA-k TODOs

## Short-term

* Documentation with full examples and quickstart, also README quickstart and remove GLPK references

## Mid-term

* Add @validate_call where possible
* Add "CompactModel" class
* Consistent argument names
* Consistent argument orders
* Add rank function
* Add in/out function in utilities, also as plot
* Add "old" bottleneck function
* Add linear-fractional programming
* Community models
* Module-based MCS (StrainDesign integration?)
* MILP gap filling
* More and better error messages
* Simple MILP EFM function for very small networks
* Switch to StringZilla

## Long-term

* Performant EFM utility
* Look up usage of populate
* numba usage?
* REPL integration (in Spyder?)
* prompt-toolkit interactive model view
* ASCII plots
* Formation energies as thermodynamic constraint alternative in Metabolite
* Cofactor swapping routines
* No-GIL support for Python ≥ 3.13
* Parallelize equilibrator handling and k_cat/k_m collection
