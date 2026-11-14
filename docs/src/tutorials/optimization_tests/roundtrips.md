# Roundtrips (`roundtrips.jl`)

## What it does

`src/analysis/roundtrips.jl` tests whether the analysis can **recover a known truth**. For one hypothesis of a setup file it injects truths for the *tested*
parameters on a grid (1 tested parameter: evenly spaced values; 2: a `round(sqrt(n))`² grid; more: a Sobol sample), generates Asimov (or toy) data at each
truth, fits it with a strategy and shows truth against fit. The grid only decides how truths are injected; the fit is free to vary every parameter that is
free in the hypothesis prior. With `--fluctuate-nuisances` (or `--toydata`) the other free parameters of the truth are redrawn from their priors, so that the
fit also has to recover fluctuated nuisance parameters.

The colour scale and the histogram show `Δχ²_mod`, either

* **`--chi2-recovery`** (default): `Δχ²_recovery = 2 (log_post_fit − log_post_truth)`. How far the fit ends from the injected truth; a diagnostic of the fit
  (a failed or incomplete optimization gives negative values), not a hypothesis-test statistic; or
* **`--chi2-mis-modeling`**: `Δχ²_mis-mod = χ²(tested parameters fixed at the truth, refit) − χ²(tested parameters free, best fit)`, the profile-likelihood
  statistic (asymptotically χ² with as many degrees of freedom as fixed tested parameters). It needs a second fit per point and roughly doubles the cost.

`--prior-effect` additionally splits `Δχ²_recovery` into `Δχ²_prior_pull = 2 (log_post_ref − truth_log_post)` (the prior moves the posterior maximum away from
the truth: expected) and `Δχ²_optimizer_residual = 2 (log_post − log_post_ref)` (≈ 0 for a well-converged fit). It uses one extra cheap local fit started at
the truth per point and is only stored, not plotted.

## How to configure it

| option | meaning |
|---|---|
| `--setup FILE` (required) | setup file with the problem |
| `--hypothesis NAME` | hypothesis to scan (default: the first one of the setup) |
| `--strategy S...` (required) | how each point is fitted: a chain of stages, e.g. `random lbfgsb` or `optim-pso:particles=20 lbfgsb:iterations=300` |
| `--nseeds N` | attempts of the whole strategy per point (the best by `log_post` is kept) |
| `--start nominal\|truth` | start point of the first stage |
| `--n-points N` | number of truths (default 25; use 100–200 for a Sobol sample) |
| `--fluctuate-nuisances`, `--toydata` | redraw the other free parameters of the truth; Poisson-fluctuated data (also shows all nuisance panels) |
| `--chi2-recovery`, `--chi2-mis-modeling`, `--prior-effect` | statistic and decomposition (see above) |
| `--objective posterior\|likelihood` | MAP or MLE fits |
| `--workers N`, `--threads T`, `--multi-node` | N worker processes with T threads each; `--multi-node` spawns them over a SLURM allocation |
| `--results-dir`, `--suffix` | output folder (default `results/<setup name>/roundtrip`) and a free tag appended to the file names |

## How to use it

```bash
julia src/analysis/roundtrips.jl --setup my_setup.jl --hypothesis theta13 \
      --strategy random lbfgsb --n-points 20 --workers 2 --suffix test

# 2D hypothesis, profile-likelihood statistic, fluctuated nuisances
julia src/analysis/roundtrips.jl --setup my_setup.jl --hypothesis theta13_dm31 \
      --strategy smseed lbfgsb --n-points 36 --chi2-mis-modeling --fluctuate-nuisances --workers 3
```

### Output

`grid_roundtrip_<asimov|toydata>_<fluctuated|nominal>_<hypothesis>_<posteriorfit|likelihoodfit>[_<suffix>]_<DDMMYYYY>.{pdf,jld2}`.

* **Figure**: for 1 tested parameter one diagonal truth-vs-fit panel; for 2 the plane of both parameters (truth and fit as separate point sets) and one
  diagonal panel per parameter; for more one diagonal panel per tested parameter. Points are coloured by `Δχ²_mod`; the marker shows whether the fit converged,
  did not converge, or has no convergence flag (global-only strategies). Failed fits (no finite result) are not drawn. With `--fluctuate-nuisances` panels for
  the nuisance parameters are added.
* **Results file** (saved before plotting, so it survives a failed plot):

```julia
using Newtrinos, FileIO
d = FileIO.load("grid_roundtrip_asimov_nominal_theta13_posteriorfit_test_14112026.jld2")
d["truth_params"], d["fit_params"]     # one NamedTuple per point (failed fits: NaN)
d["log_post"], d["truth_log_post"], d["converged"], d["n_iters"]
d["Δχ²_recovery"], d["Δχ²_mis_modeling"], d["Δχ²_prior_pull"], d["Δχ²_optimizer_residual"]   # NaN where not computed
d["n_converged"], d["n_no_convergence_flag"], d["grid_coords"], d["settings"]
```
