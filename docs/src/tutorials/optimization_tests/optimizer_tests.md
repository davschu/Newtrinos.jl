# Optimizer tests (`optimizer_tests.jl`)

## What it does

`src/analysis/optimizer_tests.jl` compares fit **strategies** (see [Optimization tests](../optimization_tests.md) for the strategy syntax and the available
optimizers) on the same problem. For one hypothesis of a setup file it injects `--n-truths` truth points (Asimov or toy data), then runs every
strategy `--nrepeats` times (each run with its own random seed) from the same blind start (the nominal parameters) and compares each run with

* the log-posterior of the **injected truth** (a successful optimizer ends at or above it; runs that end below it are counted as "below truth"), and
* a **reference fit** started *at* the truth (`--reference`, by default a long `optim-lbfgs` fit).

For every run it records the final parameters and score, the number of likelihood evaluations and the wall time (summed over *all* stages and attempts of
the strategy), the convergence flag and the trace of the best objective value so far against the number of evaluations.

## How to configure it

| option | meaning |
|---|---|
| `--setup FILE` | setup file with the problem (see [Optimization tests](../optimization_tests.md)) |
| `--hypothesis NAME` | hypothesis of the setup to test (default: the first one) |
| `--strategies S...` | strategies to compare, **one quoted argument each**, e.g. `lbfgsb "mh-de:evals=3000 lbfgsb" "optim-pso:particles=20 lbfgsb @3"` |
| `--nseeds N` | default number of attempts for strategies without their own `@N` |
| `--budget N` | default `evals` of global stages without their own (default 2000) |
| `--n-truths`, `--nrepeats` | number of injected truths (default 3) and independent runs per strategy and truth (default 1) |
| `--reference SPEC` | optimizer of the reference fit at the truth |
| `--tol X` | a run "reached the reference" if its score is within `X` of it (log-posterior units, default 0.1) |
| `--fluctuate-nuisances`, `--vary-truth-params P...` | redraw all other free parameters (or only the listed ones) of the truth from their priors |
| `--data-mode asimov\|toy` | expected values or Poisson-fluctuated data |
| `--objective posterior\|likelihood` | MAP or MLE fits |
| `--workers N` | N single-threaded worker processes (default 1: runs one after the other in this process) |
| `--no-warmup` | skip the warm-up run per optimizer and worker (without it the first run includes compilation time) |
| `--seed`, `--name`, `--results-dir` | base RNG seed, output name prefix, output folder (default `results/<setup name>/optimizer_tests`) |
| `--list-optimizers` | print the registered optimizers and their settings and exit |

## How to use it

```bash
# compare a plain local fit, two global+local chains and a repeated chain on the example setup, 3 workers
julia src/analysis/optimizer_tests.jl --setup my_setup.jl --hypothesis theta13 \
      --n-truths 3 --workers 3 --name compare \
      --strategies lbfgsb "mh-de:N=30,evals=1500 lbfgsb" "optim-pso:particles=20 lbfgsb" "random lbfgsb @5" \
      --reference "lbfgsb:iterations=200"
```

Memory: each worker needs about 2 GB for a DeepCore-sized likelihood; do not start a second multi-worker job at the same time.

### Output (file names end with the run date `DDMMYYYY`)

* `<name>_<hypothesis>_summary_<date>.txt`: one row per strategy with the number of runs, median and maximum *excess*, how many runs *reached* the
  reference, how many ended *below the truth*, median calls and time, failures and the median parameter error (relative to the prior width).
  *Excess* = reference score − score of the run: 0 means as good as the reference, positive is worse.
* `<name>_<hypothesis>_traces_<date>.pdf`: one panel per truth. y: `Δχ²_ref = 2 × excess` of the best objective value so far (symmetric log axis: linear
  inside ±0.1, where only optimizer noise is left, shown as the grey band); x: objective evaluations (log). Colour = strategy, **one line per attempt** (for `@N`
  every attempt starts again from its own beginning), line style = stage of the strategy. Dashed black line: the reference; dotted grey line: the truth's own value.
* `<name>_<hypothesis>_summary_<date>.pdf`: distance to the reference per strategy (numbered as in the legend) and wall time vs distance.
* `<name>_<hypothesis>_results_<date>.jld2`: everything, for further analysis:

```julia
using Newtrinos, FileIO
d = FileIO.load("compare_theta13_results_14112026.jld2")
d["records"]          # one NamedTuple per run: label, params, log_post, n_calls, time, converged, stages, attempt_traces, trace_x/_y/_stage, i (truth), s (strategy), ...
d["truths"], d["reference_fits"], d["excess_ref"], d["settings"], d["summary_text"]
```

Reading the results: compare strategies by *reached* (accuracy) against median *calls* and *time* (cost). A strategy with many runs *below truth* does not
find the optimum from a blind start. Global optimizers are limited by their budget and have no convergence flag; local fits that stop at their iteration limit
are reported as not converged.
