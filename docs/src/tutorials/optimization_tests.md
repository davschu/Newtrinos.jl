# Optimization tests

Fitting a global neutrino likelihood is an optimization problem with many parameters, bounded priors and expensive gradients. The scripts
in `src/analysis/` answer two questions about it:

* **Does the fit find the truth?** [Roundtrips](optimization_tests/roundtrips.md) inject a known truth (a grid or a Sobol sample of the tested
  parameters), refit and show truth against fit.
* **Which optimizer is best, and how much does it cost?** [Optimizer tests](optimization_tests/optimizer_tests.md) run several fit
  *strategies* on the same injected truths and compare accuracy, number of likelihood evaluations and wall time.

A third script, `nsi/nsi_scan.jl`, uses the same fitting engine for profile scans. All three are independent of the physics model and the data:
everything about the problem comes from a **setup file**.

## Setup files

Everything about the problem (physics, experiments, priors, hypotheses) comes from a **setup file**: a Julia file that defines `setup(; ftype)` returning
`(; name, experiments, params, priors, hypotheses)`. The scripts take it with `--setup my_setup.jl` and `--hypothesis NAME` selects one of its hypotheses.
How to write one, with a complete worked example, is described in [Creating a setup file](optimization_tests/setup_file.md); the NSI/DeepCore setup
`nsi/setups/nsi_deepcore_3y.jl` is a larger real example.

## Strategies and the available optimizers

A **strategy** is a chain of **stages**; each stage starts from the result of the previous one. On the command line a strategy is one quoted
argument `"stage stage … [@N]"`, a stage is `name` or `name:key=val,key=val`, and `@N` repeats the whole strategy N times with fresh random
draws and keeps the best (time and number of evaluations always count all attempts):

```
"lbfgsb"                                  a single local fit from the nominal start
"optim-pso:particles=50,evals=2000 lbfgsb"   global search, then a local polish
"smseed lbfgsb"                           fit the null model first, start the main fit there
"mh-de:N=30,evals=3000 lbfgsb:iterations=300 @3"   three attempts, the best is kept
```

| kind | name | what it does (settings and defaults) |
|---|---|---|
| seed | `random` | random start point drawn from the prior |
| seed | `manual:param=value,…` | start from the given parameter values |
| alias | `smseed` | local fit of the null model (tested parameters at their null values) → start of the next stage (`lbfgsb` settings) |
| alias | `mle` | `random` + local fit (`lbfgsb` settings) |
| global | `optim-pso` | particle swarm of Optim.jl (`particles=20`, `iters`, `evals`) |
| global | `optim-sa` | simulated annealing of Optim.jl (`iters`, `evals`) |
| global | `mh-pso`, `mh-de`, `mh-abc` | Metaheuristics.jl: PSO (`N=100, C1, C2, ω`), differential evolution (`N=100, F=0.7, CR=0.5`), artificial bee colony (`N=50, limit=10`); all with `evals` |
| global | `bbo`, `bbo-xnes`, `bbo-dxnes` | BlackBoxOptim.jl: adaptive DE (`method`, `population`), xNES, dxNES; all with `evals` |
| local | `lbfgsb` | L-BFGS-B through BAT's `OptimizationAlg` (`iterations=1000, reltol=1e-7, m=10, ad`); the recommended local stage |
| local | `optim-lbfgs`, `optim-gradientdescent`, `optim-conjugategradient`, `optim-neldermead` | Optim.jl local methods (`iterations=2000, g_tol=1e-6, f_tol, x_tol, step_size, ad`) |
| local | `newtrinos` | `Newtrinos.find_mle` (no settings; full fit) |
| local | `ext` | `Newtrinos.find_mle_ext` (`iterations`, `g_tol`) |

Rules and remarks:

* Seed, alias and global stages can only be the **first** stage; local stages can follow anywhere.
* Global stages stop after `evals` objective evaluations (default: the script's `--budget`), so they have no convergence flag (`converged = missing`).
  A strategy that is only a global search is not polished.
* From a blind start plain `optim-lbfgs` can jump onto the prior bounds and stall; use `lbfgsb`.
* The start point of every strategy is the nominal parameters of the setup, not the truth.
* `ad` selects the automatic-differentiation backend (`auto`, `forward`, `polyester`); inside multi-threaded runs `forward` is used.
* `julia src/analysis/optimizer_tests.jl --list-optimizers` prints the registry with all settings.
* Own optimizers can be added with `register_optimizer!(name, kind, fn, defaults)` in `optimizer_tests_common.jl`.

## Using the engine from Julia

The same engine can be used outside the scripts:

```julia
include("src/analysis/optimizer_tests_common.jl")

# one fit with a strategy (returns the record of the last stage plus the cost and traces of the whole strategy)
r = find_mle_advanced(likelihood, prior, params; strategy = ["mh-de:evals=2000", "lbfgsb"], nseeds = 3, null_prior = null_prior)
r.params, r.log_post, r.total_calls, r.total_time, r.stages, r.attempt_traces

# a profile scan with a strategy, optionally warm-started along the grid; returns a Newtrinos.NewtrinosResult
res = profile_advanced(likelihood, priors, (; θ₁₃ = 20), params; strategy = ["random", "lbfgsb"], sequential = true)
res.meta["optimizer_runs"]   # converged, n_iterations, time, n_calls, n_attempts, stages, trace per grid point
```

## Output files

All files written by `optimizer_tests.jl`, `roundtrips.jl` and `nsi_scan.jl` end with the date of the run (`DDMMYYYY`, e.g. `14112026`) before the
extension, so repeated runs do not overwrite each other. Results are saved with `FileIO.save` as `.jld2` (load with `using Newtrinos, FileIO`);
profile scans are `Newtrinos.NewtrinosResult` objects under the key `"result"`.

## Parallelization

* `--workers N` starts N single-threaded worker processes (`optimizer_tests.jl`, `roundtrips.jl`, `nsi_scan.jl`); timings stay comparable. Each worker
  needs about 2 GB of memory for a DeepCore-sized likelihood. Run one multi-worker job at a time.
* `--threads T` (`roundtrips.jl`, `nsi_scan.jl`) sets the threads per worker; with `--workers 1` the grid runs in the calling process and uses the
  threads of `julia -t T`.
* Results do not depend on the number of workers or threads: every grid point, truth and run has its own random stream.
* `--multi-node` (`roundtrips.jl`, `nsi_scan.jl`) spawns the workers across a SLURM allocation.
