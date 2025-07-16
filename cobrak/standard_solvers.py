from platform import platform

from joblib import cpu_count

from .dataclasses import Solver

BARON = Solver(
    name="gams",
    solver_options={
        "solver": "baron",
    },
    solve_extra_options={
        "add_options": [
            "GAMS_MODEL.optfile = 1;",
            "$onecho > baron.opt",
            "$offecho",
        ]
    },
)

CPLEX = Solver(
    name="cplex_direct",
)

CPLEX_FOR_VARIABILITY_ANALYSIS = Solver(
    name="cplex_direct",
    solver_options={
        "threads": 1,
        "lpmethod": 1,
    },
)

GLPK = Solver(
    name="glpk",
    solver_attrs={"_version_timeout": 180},
)

GUROBI = Solver(
    name="gurobi_direct",
)

GUROBI_FOR_VARIABILITY_ANALYSIS = Solver(
    name="gurobi_direct",
    solver_options={
        "threads": 1,
    },
)

hsllib = "/scratch/CoinHSL/coinhsl-2024.05.15/BUILD/lib/x86_64-linux-gnu/libcoinhsl.so"
if platform().startswith("macOS"):
    hsllib = "/Users/pbekiaris/0Programs/CoinHSL.v2024.5.15.aarch64-apple-darwin-libgfortran5/lib/libcoinhsl.dylib"
elif cpu_count() > 16:
    hsllib = "/u/pbekiaris/CoinHSL_2024_BUILT/lib/x86_64-linux-gnu/libcoinhsl.so"
elif cpu_count() > 6:
    hsllib = (
        "/mechthild/home/bekiaris/CoinHSL_2024_BUILT/lib/x86_64-linux-gnu/libcoinhsl.so"
    )

IPOPT_MA57 = Solver(
    name="ipopt",
    solver_options={
        "max_iter": 4_000,
        "halt_on_ampl_error": "yes",
        "mu_strategy": "adaptive",
        "corrector_type": "primal-dual",
        "acceptable_tol": 1e-8,  # same as default "tol"
        "acceptable_constr_viol_tol": 0.0001,  # same as default "constr_viol_tol"
        "acceptable_dual_inf_tol": 1.0,  # same as default "dual_inf_tol"
        "acceptable_iter": 0,
        "linear_solver": "ma57",
        "hsllib": hsllib,
        "nlp_scaling_method": "none",
        "ma57_automatic_scaling": "no",
        "ma57_block_size": "32",
        "ma97_small": 1e-30,
        "ma97_u": 1e-6,
    },
)

IPOPT = Solver(
    name="ipopt",
    solver_options={
        "max_iter": 4_000,
        "halt_on_ampl_error": "yes",
    },
)

IPOPT_LONGRUN = Solver(
    name="ipopt",
    solver_options={
        "max_iter": 100_000,
        "halt_on_ampl_error": "yes",
        "mu_strategy": "adaptive",
        "corrector_type": "primal-dual",
        "linear_solver": "ma57",
        "hsllib": hsllib,
        "acceptable_tol": 1e-8,  # same as default "tol"
        "acceptable_constr_viol_tol": 0.0001,  # same as default "constr_viol_tol"
        "acceptable_dual_inf_tol": 1.0,  # same as default "dual_inf_tol"
        "acceptable_iter": 0,
        "nlp_scaling_method": "none",
        "ma57_automatic_scaling": "no",
    },
)

SCIP = Solver(
    name="scip",
    solver_attrs={"_version_timeout": 180},
)
