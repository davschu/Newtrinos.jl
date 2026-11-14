# nsi/: profile and scan driver and the NSI setup

| file | what it is |
|---|---|
| `nsi_scan.jl` | profile (or plain likelihood scan) of one hypothesis of a setup file, fitted with a strategy at every grid point |
| `setups/nsi_deepcore_3y.jl` | NSI (Standard parametrization) with the DeepCore 3-year high-statistics sample: hypotheses `eμ`, `eτ`, `μτ`, `ee_μμ`, `ττ_μμ`, `all` |
| `deepcore_3year_analysis/` | the frozen scripts and results of the first 3-year analysis (git tag `nsi-3y-scripts`) and old SLURM job scripts; not used by the current scripts |

The engine (strategies, optimizers, `profile_advanced`) and the setup-file format are described in
[`docs/src/README_optimizer_tests.md`](../docs/src/README_optimizer_tests.md). `nsi_scan.jl` itself is not NSI-specific: it works with any `--setup` file
(see the docs page *Creating a setup file*).

## nsi_scan.jl

```bash
julia nsi/nsi_scan.jl --setup nsi/setups/nsi_deepcore_3y.jl --name run1 --task profile --hypothesis ττ_μμ \
      --objective posterior --grid-points 25 --strategy "optim-pso:particles=50,iters=30" lbfgsb --workers 4
```

* `--task profile`: the hypothesis' tested parameters are scanned and every other free parameter is fitted at each grid point
  (`Newtrinos.profile` with a strategy). `--task scan`: only the likelihood at the grid points, nothing is fitted.
* Grid: 1 tested parameter → `--grid-points` evenly spaced values; 2 → `round(sqrt(n))` per dimension; more than 2 → `n` Sobol points.
* `--strategy STAGE...` (required): the fit at every grid point, a chain of stages (see the strategy table in `docs/src/README_optimizer_tests.md`).
  A **seed or alias first stage** (`random`, `manual`, `smseed`, `mle`) does not depend on the grid point; it is computed **once** and shared
  by all points (the scanned parameters are overwritten per point). A **global first stage** runs at every grid point.
* `--nseeds K`: K independent attempts per grid point (for a seed first stage: K different seeds); the best by log-posterior is kept.
  `--seed` is the base RNG seed (every grid point has its own stream, so results do not depend on the number of workers or threads).
* `--objective posterior|likelihood`: MAP or MLE fits.
* `--workers N` (single-threaded workers by default, `--threads` per worker): grid points are distributed with `pmap`; with `--workers 1` the grid
  runs in the calling process, threaded if Julia was started with `-t`. `--multi-node` spawns the workers across a SLURM allocation
  (`--workers` must equal `SLURM_NTASKS`).
* `--results-dir` (default `results/<setup name>/scan`), `--name`.

### Output

One file per run, **`<name>_profile_<hypothesis>_<n>points_<stages>_nseeds<k>_<DDMMYYYY>.jld2`** (task `scan`:
`<name>_scan_<hypothesis>_<n>points_<DDMMYYYY>.jld2`), holding a `Newtrinos.NewtrinosResult` under the key `"result"`:

```julia
using Newtrinos, FileIO
r = FileIO.load("run1_profile_ττ_μμ_25points_optim-pso+lbfgsb_nseeds1_14112026.jld2")["result"]
r.axes            # grid axes, e.g. (Δ_τμ = [...],); Sobol: (point = 1:N,)
r.values.llh, r.values.log_posterior, r.values.θ₂₃, ...   # one array per parameter (shaped like the grid)
r.meta["optimizer_runs"]   # converged, n_iterations, time, n_calls, n_attempts, stages, trace_x/trace_y/trace_stage, strategy
r.meta                     # hypothesis, setup, objective, strategy, nseeds, seed, exec_time, scan_keys, ... (Sobol: "sobol" => true)
```

`meta["optimizer_runs"]` arrays have the shape of the grid; `time` and `n_calls` count all stages and attempts of the point. The `.jld2` can
only be loaded with `Newtrinos` loaded (it stores the `NewtrinosResult` type).
