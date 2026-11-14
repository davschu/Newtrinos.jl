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

## Helper functions

All helpers live in `src/analysis/optimizer_tests_common.jl`, which the scripts `include` (on the workers too). Functions starting with `_` are internal.
The ones you may call yourself, grouped by purpose:

**Setup and hypotheses**

| function | purpose |
|---|---|
| `finalize_setup(S)` | checks the setup returned by `setup()` (fields, every hypothesis, tested parameters free) and returns it |
| `hypothesis_info(S, name)` | everything derived from a hypothesis: prior, null values, tested parameters, null prior, nuisance parameters |

**Injected truths**

| function | purpose |
|---|---|
| `truth_grid(prior, params, tested; n_points)` | truths over the tested parameters: 1D grid, 2D grid or Sobol sample, depending on how many are tested |
| `controlled_truth(rng, params, prior, tested)` | one random truth: only the tested parameters are drawn from the prior, the rest stays nominal; fixed parameters are forced to their fixed value |
| `fluctuate_nuisances(rng, truth, prior, names)` | redraws the nuisance parameters of a truth from their priors (`--fluctuate-nuisances`) |
| `randomize_free_params(rng, params, prior)` | draws all free parameters from the prior and leaves fixed ones alone; also works for parameters fixed by `Newtrinos.condition` |

**Priors and objective**

| function | purpose |
|---|---|
| `extract_bounds(prior, params)` | splits a prior into fixed and free parameters and returns the bounds and start vector of the free ones |
| `flat_prior(prior)` | replaces every free prior by a uniform one over its range (for the `likelihood` objective) |
| `make_objective(likelihood, prior, params; objective)` | the function to minimize: negative log posterior or log likelihood over the free parameters |
| `CountingObjective(f, budget)` | wraps an objective, counts evaluations, stops at the budget and records the best value so far (the traces) |

**Optimizers and strategies**

| function | purpose |
|---|---|
| `register_optimizer!(name, kind, fn, defaults)` | adds an optimizer to the registry `OPTIMIZERS` |
| `parse_optimizer_spec("name:key=val")` | splits a stage into name and settings (typed from the defaults) |
| `parse_strategy_arg("stage stage @N")` | turns one command-line strategy argument into stages and number of attempts |
| `normalize_strategy(strategy)` | checks a strategy (stage order, known names) and expands aliases such as `smseed` and `mle` |
| `run_optimizer(spec, likelihood, prior, params; budget, rng)` | runs one registered optimizer (a single stage) |
| `global_find_mle(...)` / `local_find_mle(...)` | the runners behind global and local stages; `local_find_mle` evaluates the point if no parameter is free |
| `find_mle_advanced(likelihood, prior, params; strategy, nseeds, ...)` | the main fit: runs a whole strategy (all attempts) and returns the best result with total calls, time and traces |

**Profiles**

| function | purpose |
|---|---|
| `profile_advanced(likelihood, priors, vars, params; strategy, ...)` | `Newtrinos.profile` with a strategy at every grid point; returns a `NewtrinosResult` |
| `_profile_plan`, `_profile_point_fit`, `_profile_result` | its three steps (plan the points, fit one point, assemble the result); `nsi/nsi_scan.jl` calls them directly to run the points with `pmap` |

For which script calls which function, see `docs/src/README_optimizer_tests.md`.

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
