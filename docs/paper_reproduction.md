# Reproduce publication results

!!! warning "Solver usage"
    The iCH360 scripts mentioned herein usually use CPLEX as LP solver. Please refer to this documentation's installation chapter for instructions on how to install it for COBRA-k. Alternatively, simply switch all instances where "CPLEX" is called to another solver (e.g. the pre-installed SCIP) or add a ```solver=YOUR_SOLVER``` option to any case where an optimization is called.

    Also, for the iCH360 calculations, the IPOPT sub-solver MA57 is used, such that the ```IPOPT_MA57``` solver with personal HSLLIB path settings from the publication's authors is loaded from ```cobrak.standard_solvers```. To use MA57, obtain it through HSLLIB from here (free for academics):

    [https://licences.stfc.ac.uk/product/coin-hsl](https://licences.stfc.ac.uk/product/coin-hsl)

    ...and change the ```"hsllib"``` ```solver_option```in ```IPOPT_MA57``` to your HSLLIB path.

!!! warning "Re-run of existing calculation"
    Many iCH360 calculations will not run because existing calculation results are found (this prevents unneccessary double calculations). To mitigate this problem, simply delete all folders in the ```iCH360``` subfolder except for the ```external_resources``` folder.

    *Note:* For the $k_{cat}$ and $K_M$ variation analysis, you have to un-zip the three zip files there ("used_cobrak_models.zip", "variability_dicts.zip" and "best_evolution_results.zip") into the variation results folder beforehand if you want to get exactly the same random variation calculations as the ones shown in COBRA-k's publication.

To reproduce the results of COBRA-k's initial publication:

* Bekiaris & Klamt, COBRA-k: A powerful framework bridging constraint-based and kinetic metabolic modeling. <i>Sci. Adv.</i> 12, eaeb3022 (2026). [DOI:10.1126/sciadv.aeb3022](https://doi.org/10.1126/sciadv.aeb3022)

...head over to the ```examples``` subfolder in COBRA-k's main folder. There, you can already find the first important Python scripts that pre-processed some data:

* FIRST_A_read_out_sabio_rk.py: Reads out SABIO-RK and creates a cache into the ```common_needed_external_resources```folder
* FIRST_B_get_ec_number_transfers.py: Reads out Expasy obsoleted EC number data and creates a cache into the ```common_needed_external_resources``` folder
* FIRST_C_read_out_met_concs_file.py: Reads out BiGG metabolites file obsoleted EC number data and puts a more machine-readable version into the ```common_needed_external_resources``` folder
* FIRST_D_read_out_enzyme_abundance_file.py: (Specific for the publication) Reads out the *in vivo* enzyme data

So far, all these steps were just some caching steps to make our subsequent model creation steps reproducibly faster.

For the toy model calculations, simply run ```examples/toymodel/run_toymodel_calculations.py```.

Let's continue with ```iCH360``` (regarding the iCH360 COBRA-k variant described in COBRA-k's publication) in the appropriate subfolder. There, simply follow the alphabetically sorted scripts and run them with ```uv run```. The only more complex script for usage is ```H_run_calculations.py```. As input, it takes JSON with a reproducible run configuration (with model settings and so on). To generate the run setting JSONs of COBRA-k's publication on your local machine, run ```main_paper_calculations.py local``` in COBRA-k's main folder. Note the ```local``` option, which creates the JSONs in a ```main_paper_calculations_jsons``` subfolder. Without the ```local``` option, you can use to run all the main paper calculations on a SLURM computer cluster as it creates a SLURM file and ```sbatch```es it for each JSON run config.

If you just want to play around with iCH360_cobrak, you can simply import it through ```from cobrak.example_models import iCH360_cobrak```. iCH360_cobrak is now a ```Model``` variable with the full model.
