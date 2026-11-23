# Analysis scripts: optimizer comparison and roundtrips

Physics- and dataset-independent tooling to **test how well fits recover known truths** and **compare optimizers**. Everything
about the problem (physics, data, priors, hypotheses) comes from a *setup file*; the scripts contain no model-specific code.
The user-facing documentation with more examples is the docs page *Tutorials → Optimization tests*.

| file | what it is |
|---|---|
| `optimizer_tests.jl` | compares optimizer **strategies** on injected truths (accuracy, cost, convergence traces) |
| `roundtrips.jl` | inject truths on a grid / Sobol sample, refit, and plot truth vs fit (`Δχ²_recovery` or `Δχ²_mis-mod`) |
| `optimizer_tests_common.jl` | the engine both scripts (and `nsi/nsi_scan.jl`) share: optimizer registry, strategies, `find_mle_advanced`, `profile_advanced`, truth/setup helpers |
| `nsi/setups/nsi_deepcore_3y.jl` | the NSI / DeepCore 3-year setup (see `nsi/README.md`) |

Run from the repository root (`my_setup.jl` is a setup file, see below and the docs page *Creating a setup file*), e.g.

```bash
julia src/analysis/optimizer_tests.jl --setup my_setup.jl --hypothesis theta13 \
      --strategies lbfgsb "mh-de lbfgsb" "random lbfgsb @3" --workers 3
julia src/analysis/roundtrips.jl --setup my_setup.jl --hypothesis theta13 \
      --strategy random lbfgsb --n-points 10 --workers 2
```

## Setup files

A setup file defines `setup(; ftype)` returning `(; name, experiments, params, priors, hypotheses)`; `hypotheses` is an ordered
dictionary `name => (; prior, null)`.

* `experiments`: NamedTuple of configured experiments; `params`/`priors`: from `Newtrinos.get_params/get_priors(experiments)`.
* A **hypothesis** is one question "is effect X present and can the fit recover it?":
  * `prior`: the prior in which the parameters of X are **free** and every parameter that does not belong to the analysis is
    **fixed**. Build it with `Newtrinos.condition(priors, Dict(:name => value, ...), params)` (listed parameters become constants).
  * `null`: NamedTuple of the **tested** parameters with their values when X is absent, e.g. `(ε_eμ_abs = 0.0, δ_eμ = 0.0)`.
* From these the scripts derive: tested parameters = `keys(null)` (their order is the order of the scan axes); the null model
  (`condition(prior, null, params)`, used by the `smseed` stage); the scan type (1 tested parameter → 1D grid, 2 → 2D grid, more →
  Sobol sample); the nuisance parameters = all other free parameters of `prior` (fitted in every roundtrip, redrawn by
  `--fluctuate-nuisances`).
* Fix a nuisance parameter: `priors = Newtrinos.condition(priors, Dict(:atm_flux_updown_sigma => 0.0), params)`. Change a prior:
  `@reset priors.deepcore_lifetime = Uniform(0, 3.8)`. Do both on `priors` **before** building the hypotheses.
* Every tested parameter must be free in the hypothesis prior (checked by `finalize_setup`). If *all* free parameters are tested ones,
  a profile fit has nothing left to fit; the point is then just evaluated.

## Strategies and optimizers

A **strategy** is a chain of **stages**, each started from the result of the previous one. In the command line a strategy is one quoted
argument `"stage stage ... [@N]"`; a stage is `name` or `name:key=val,key=val`; `@N` runs the whole strategy N times with fresh random
draws and keeps the best (time and calls always count all attempts). `optimizer_tests.jl --list-optimizers` prints the registry.

| kind | stages | rule |
|---|---|---|
| seed | `random` (random start from the prior), `manual:param=value,...` | first stage only |
| alias | `smseed` (local fit of the null model → start), `mle` (random + local fit) | first stage only |
| global | `optim-pso`, `optim-sa`, `mh-pso`, `mh-de`, `mh-abc`, `bbo`, `bbo-xnes`, `bbo-dxnes` | first stage only; stop at `evals` objective evaluations |
| local | `lbfgsb`, `optim-lbfgs`, `optim-gradientdescent`, `optim-conjugategradient`, `optim-neldermead`, `newtrinos`, `ext` | anywhere after the first |

Notes: a strategy that is a single global search is not polished and has no convergence flag (`converged = missing`). Plain
`optim-lbfgs` can collapse onto the prior bounds from a blind start; use `lbfgsb` (L-BFGS-B, tunable: `iterations`, `reltol`, `m`, `ad`).
`newtrinos` and `ext` do not take `iterations` and run full fits. The start point of every blind strategy is the setup's nominal
parameters (not the truth).

## optimizer_tests.jl

Compares strategies (`--strategies`, default `optim-lbfgs`, `optim-pso:particles=20`, `mh-de`). For each of `--n-truths` injected truth
points of one hypothesis (data: `--data-mode asimov|toy`), each strategy is run `--nrepeats` times (own seed each) and compared with

* the log-posterior of the injected truth (a successful optimizer ends at or above it; "below truth" counts failures), and
* a truth-seeded reference fit (`--reference`, started at the truth).

Recorded per run: parameters, score, `n_calls` and `time` (summed over all stages and attempts), the trace of the best objective so far
against evaluations (per attempt and combined), converged flag. With `--workers N` the runs are distributed over N single-threaded
worker processes (timings stay comparable; every worker needs about 2 GB). A warm-up run per optimizer and worker removes
compilation time from the timings (`--no-warmup` disables it).

Options: `--setup`, `--hypothesis` (default: the setup's first), `--strategies`, `--nseeds` (default attempts for strategies without `@N`),
`--budget` (default `evals` of global stages), `--n-truths`, `--nrepeats`, `--reference`, `--tol` (a run "reached" the reference if its
excess is below it), `--fluctuate-nuisances`, `--vary-truth-params`, `--data-mode`, `--objective posterior|likelihood`, `--seed`, `--name`,
`--results-dir` (default `results/<setup name>/optimizer_tests`), `--no-warmup`, `--workers`, `--list-optimizers`.

Output (every file name ends with the run date `DDMMYYYY`):

| file | content |
|---|---|
| `<name>_<hypothesis>_summary_<date>.txt` | table per strategy: runs, median/max excess, reached, below truth, median calls and time, failures, parameter error |
| `<name>_<hypothesis>_results_<date>.jld2` | `records` (every run), `truths`, `reference_fits`, `excess_ref`, `excess_truth`, `settings`, `summary` |
| `<name>_<hypothesis>_traces_<date>.pdf` | `Δχ²_ref` of the best objective so far vs objective evaluations, one panel per truth; colour = strategy, one line per attempt, line style = stage; grey band = optimizer noise (`|Δχ²| ≤ 0.1`) |
| `<name>_<hypothesis>_summary_<date>.pdf` | distance to the reference per strategy, and wall time vs distance (only written when more than one strategy is compared) |

`excess` = reference score − score of the run (0: as good as the reference; > 0 worse). The plots show `Δχ²_ref = 2 × excess` on a
symmetric log axis (linear inside ±0.1).

## roundtrips.jl

Injects truths, refits and shows how well the truth is recovered. The **grid** controls only how truths are injected (via `truth_grid`):
1 tested parameter → evenly spaced values, 2 → `round(sqrt(n))` per dimension, more → Sobol sample (use `--n-points` 100–200). The fit is
still free to vary every non-fixed parameter of the hypothesis prior. With `--fluctuate-nuisances` (or `--toydata`) the other free
parameters of the truth are redrawn from their priors, and per-parameter panels are added.

Per hypothesis the figure shows (by the number of tested parameters): 1 → one truth-vs-fit diagonal panel; 2 → the plane of both
parameters (truth and fit as separate point sets) plus one diagonal panel per parameter; more → one diagonal panel per tested parameter.
Colour/histogram: `Δχ²_mod`, one of

* `--chi2-recovery` (default): `Δχ²_recovery = 2 (log_post_fit − log_post_truth)`, i.e. how far the fit ends from the injected truth.
  A pipeline/optimizer diagnostic ("did the fit find the truth"), not a hypothesis-test statistic.
* `--chi2-mis-modeling`: `Δχ²_mis-mod = χ²(tested parameters fixed at the truth, refit) − χ²(tested parameters free, best fit)`: the profile-likelihood
  statistic, asymptotically χ²-distributed with dof = number of fixed tested parameters (Wilks). Needs a second fit per point (≈ doubles the cost).

`--prior-effect` decomposes `Δχ²_recovery` into `Δχ²_prior_pull = 2 (log_post_ref − truth_log_post)` (the prior can pull the posterior maximum
away from the truth; expected, not a bug) and `Δχ²_optimizer_residual = 2 (log_post − log_post_ref)` (should be ≈ 0 for a well converged fit),
using one extra cheap local fit started at the truth. Both are NaN without the flag, and neither is drawn.

Options: `--setup` (required), `--hypothesis`, `--results-dir` (default `results/<setup name>/roundtrip`), `--n-points`, `--strategy`
(required, as above), `--nseeds`, `--start nominal|truth`, `--objective`, `--fluctuate-nuisances`, `--toydata`, `--chi2-recovery`,
`--chi2-mis-modeling`, `--prior-effect`, `--suffix`, `--workers`, `--threads` (per worker), `--multi-node` (SLURM).

Output: `grid_roundtrip_<asimov|toydata>_<fluctuated|nominal>_<hypothesis>_<posteriorfit|likelihoodfit>[_<suffix>]_<date>.{pdf,jld2}`. The
`.jld2` is saved **before** plotting and holds `grid_coords`, `truth_params`, `fit_params`, `log_post`, `truth_log_post`, `converged`, `n_iters`,
the `Δχ²_*` arrays, timings and the `settings`. A failed fit has `NaN` parameters (they are not drawn); how many fits failed can be read
from the results file only. Failed fits that are drawn would be red crosses.

## sensitivity.jl

Sensitivity studies on one hypothesis of a setup file. The truth is the setup's parameters (pinned parameters of the hypothesis forced);
the parameters studied are the free parameters of the hypothesis prior, the tested parameters are the parameters of interest. `--method`:

* `fisher` (no fits): the Hessian of the log posterior at the truth. `--information expected` (default) uses the Fisher information of the
  Poisson model, the Hessian of `Σ μ₀ log μ(q) − μ(q)` with the expected counts μ₀ at the truth fixed (`J' diag(1/μ) J`, no data, no rounding,
  positive semi-definite); `--information observed` the Hessian of the log likelihood of the rounded Asimov counts, which can be indefinite
  where the data are insensitive (such truths get no result). From it: the expected
  uncertainty of every parameter with all others free (`sigma_marg`) and fixed (`sigma_cond`) against its prior width (which parameters the
  data constrain), the impact of every nuisance on the tested parameters (first-order shift for a +1σ move, and `sigma` with that nuisance
  fixed), and the **a priori correlation matrix** (`Newtrinos.correlation` at the truth). Where the data are injected is chosen with `--truth`:
  `nominal` (default, the setup's parameters), `grid` (`--n-truths N` points of `truth_grid` over the whole prior ranges of the tested parameters, including the bounds:
  global sensitivity; 1 tested parameter: a line, 2: a grid, more: Sobol points), `range:lo,hi` (one tested parameter) or `range:name=lo,hi;name2=lo,hi` (the same
  grid restricted to intervals, e.g. a suspicious region: `--truth range:-1.3,-0.7 --n-truths 13`; tested parameters not listed keep their whole range),
  `point:name=x,name2=y` (ONE truth: these parameters, tested or nuisances, set on top of the nominal ones) and `bf:result.jld2` (ONE truth: the best fit of a
  `NewtrinosResult` stored under the key `"result"`, as written by `save_result`; parameters the result does not contain stay nominal, with a warning). The intervals and values must lie in the prior
  range. The summary states which truth was used (an a priori forecast at the nominal truth, or post-fit at a best fit) and the `.jld2` keeps it in
  `settings["truth"]`. With `--drop-boundaries` the grid is cell-centred (`lo + (i-1/2)/n * width`), so no truth lies on a bound or at a zero-effect point such as a
  phase with zero magnitude. A Gaussian approximation: use `scan` where the likelihood is strongly non-quadratic or multimodal (δCP, θ₂₃ octant), and be careful with truths where
  the model changes quickly (a vacuum-like truth such as Δ_eμ = −1 for NSI: the matter effect switches on there, sigma is very small and correlations flip sign). If the matrix is singular at the
  nominal truth (e.g. a phase whose magnitude is zero), use `--truth grid`.
* `scan`: one-at-a-time scans of every free scalar parameter (or `--params`) over its **full prior range** (`--n-points` points plus the truth itself) on Asimov data:
  with all other parameters fixed (likelihood evaluations only) and, unless `--no-profile`, with all other free parameters re-fitted
  (`--strategy`) at every point. `Δχ² = 2 (llh(truth) − llh(point))`; the summary lists the Δχ² at the worst grid point, the 1σ interval
  (Δχ² = 1) and, for the profiled scans, the largest shift of the tested parameters (the impact of a nuisance on the result).
* `toys`: `--n-toys` toy data sets injected at the **null** of the hypothesis; per toy a null fit and then a full fit started from the null fit (so q ≥ 0 up to the optimizer tolerance; use local strategies, a seed or global first stage would ignore the start); the distribution of
  `q = 2 (log_post_full − log_post_null)` (without the tested parameters' own prior density, which differs between the two fits, e.g. −ln(width) of a Uniform) is compared with `Chisq(number of tested parameters)` (Wilks check, quantile table and histogram).

Options: `--setup` (required), `--hypothesis`, `--method`, `--params`, `--n-points`, `--no-profile`, `--strategy` (default `lbfgsb`; ONE strategy, its stages as separate arguments or in one quoted string, e.g. `"mle optim-lbfgs:iterations=2000"`;
unlike `roundtrips.jl` several arguments are not several strategies; a last token `@N`, `scan` only, runs N independent attempts per fit and keeps the best),
`--n-truths`, `--truth`, `--information`, `--drop-boundaries`, `--top-n`, `--n-toys`, `--seed`, `--name`, `--results-dir` (default `results/<setup name>/sensitivity`), `--workers`.
Output: `<name>_<hypothesis>_<method>_{summary,results}_<date>.{txt,jld2}` plus the plot. For `fisher` one figure, `_overview`
(PDF), with the truths resolved instead of averaged: row 1 σ / prior width of every tested parameter against its own truth (σ_marg solid, σ_cond dashed, log axis, a dotted line at
σ = width, where the data no longer constrain the parameter within its prior range; several tested parameters: one panel each, markers) next to the nuisances (σ_marg and σ_cond over
the prior width, mean over the truths with min-max whiskers); row 2 the RMS of the correlation matrix over the truths (sign not shown, a smell test) next to
the range over the truths of **impact / σ_marg** of each nuisance (shift of the tested parameter for a +1σ_marg move of the nuisance, in units of the σ_marg of the tested parameter,
equal to their correlation; a bar crossing zero is a sign flip, which the RMS matrix cannot show). At most `--top-n` (default 25) nuisances are drawn, the ones with the largest |impact| at any
truth (the title says so); the summary text and the `.jld2` hold all of them, and the height of the panels follows the number of drawn rows. With ONE truth (`point:`, `bf:`, `nominal`)
there is nothing to resolve over truths, so the overview shows the actual numbers instead: the σ / prior width of the tested parameter(s) (listed first) and the nuisances in one bar panel, the signed
correlation matrix at the truth and the signed impact / σ_marg bars. Truths without a result are named in the title; one figure for `scan` and `toys`.
The `.jld2` and the summary are written after the plot, also if the plot fails (the error is logged). Plot code: `src/analysis/sensitivity_plots.jl`
(`sensitivity_data(results)`, `static_overview(D; top_n)`), usable on a loaded `.jld2`.

**Interactive explorer (on demand).** `nsi/deepcore_8year_analysis/results/sensitivity_results/make_explorer_html.jl` turns the results `.jld2` of a `fisher` run into ONE self-contained
HTML file with a slider (and the arrow keys) over the injected truth; it is not written by `sensitivity.jl` (the file is large):
`julia --project=<Newtrinos.jl> make_explorer_html.jl <results.jld2> [--out explorer.html] [--top-n 25] [--max-frames 60]`. Every truth is a pre-rendered frame (about 0.4 MB; at most `--max-frames`,
evenly subsampled) with the σ / prior width of the tested parameter(s) against the truth (the selected one marked) and of the nuisances at the truth, and below the correlation matrix and the
impact / σ_marg of the nuisances at the truth. No Julia server or kernel is needed to open it; `page.html#N` starts at the (N+1)-th frame.

The functions behind it work on any result: `Newtrinos.correlation(result, likelihood)` (correlation at a profile's best fit) and
`Newtrinos.gof_pvalue(result, experiments)` (Poisson goodness of fit per experiment) are documented in the API page "Analysis API".

## Who calls what

```
optimizer_tests.jl ─┐
roundtrips.jl ──────┼─ finalize_setup, hypothesis_info                 setup file → problem, tested/null/nuisance parameters
nsi/nsi_scan.jl ────┘
optimizer_tests.jl, roundtrips.jl ── truth_grid / controlled_truth / fluctuate_nuisances   (truth injection)

optimizer_tests.jl ─ parse_strategy_arg ─┐
roundtrips.jl, nsi_scan.jl ─ normalize_strategy ─┴─► find_mle_advanced(strategy, nseeds, starts)
                                                         └─ _run_strategy ─ per stage ─ run_optimizer ─ OPTIMIZERS[name].fn
                                                              seed/alias:  _seed_random, _seed_manual, _expand_smseed, _expand_mle
                                                              global:      global_find_mle ─ CountingObjective ─ _backend_{optim_pso, optim_sa, mh_*, bbo}
                                                              local:       local_find_mle (optim | lbfgsb | newtrinos | ext)
                                                         └─ traces: _stage_trace, _strategy_trace (per attempt and combined)

profile_advanced ─ _profile_plan ─ _profile_point_fit ─ find_mle_advanced      (also called piecewise by nsi_scan.jl)
                 └ _profile_result ─ _optimizer_runs ─ Newtrinos.assemble_profile_results → NewtrinosResult
                 └ _warm_params_from_results, Newtrinos._sequential_order       (warm / sequential starts)
roundtrips.jl ─ local_find_mle (extra fits: --prior-effect, --chi2-mis-modeling), with_local_defaults, final_local_spec
```

`find_mle_advanced` is also the entry point for use in your own scripts:

```julia
include("src/analysis/optimizer_tests_common.jl")
r = find_mle_advanced(likelihood, prior, params; strategy = ["mh-de:evals=2000", "lbfgsb"], nseeds = 3, null_prior = null_prior)
r.params, r.log_post, r.total_calls, r.stages, r.attempt_traces
res = profile_advanced(likelihood, priors, (; θ₁₃ = 20), params; strategy = ["random", "lbfgsb"], sequential = true)   # NewtrinosResult
res.meta["optimizer_runs"]    # per grid point: converged, n_iterations, time, n_calls, n_attempts, stages, trace_x/_y/_stage
```

## Implementation notes

Things that look odd in the code and why they are there:

* **Seeds:** every grid point, truth and attempt has its own random stream (`Xoshiro(seed + i)` in profiles, `Xoshiro(1000 k + s)` in roundtrips), so results do
  not depend on the number of workers or threads. `Random.seed!(rng)` without an integer is *not* reproducible; toy data (`generate_toy_data`) draws from the
  global RNG, which is therefore seeded with an integer per grid point.
* **Start point:** optimizers must start from `_fix_conditioned(p, prior)`, not from the raw nominal parameters: `bat_findmode`'s `ExplicitInit` rejects a start
  that differs from the constants a conditioned prior holds fixed (`ArgumentError: Cannot set constant value to a different value`). The same applies when the
  scanned parameters of a seed are overwritten per grid point.
* **Prior bounds:** a start exactly on a bound makes BAT's transformation non-finite, so local fits nudge the start `1e-8` inside the bounds (`_nudge_into_bounds`).
  Note that `smseed` leaves tested parameters `1e-8` from a bound, so the main fit can stall there.
* **Workers:** a plain `@everywhere function` does not see driver-side variables; shared data is bundled in one NamedTuple (`ctx`, roundtrips) or broadcast
  once (`@everywhere plan = $plan`, `nsi_scan.jl`). In `nsi_scan.jl` the likelihood is a worker-global used inside `pmap` instead of a call to `profile_advanced`,
  because a closure in a library function would serialize the whole likelihood to every worker. `using Accessors` must be loaded in a separate `@everywhere`
  block before the file that uses `@reset` (world age).
* **Threads vs processes:** Julia's garbage collector stops all threads of a process, so a few multi-threaded worker processes beat one process with many threads.
* **AD backend:** with threads/workers the local fits use `ad = :forward` (PolyesterForwardDiff inside already-parallel runs is slow and memory hungry).
* **Plots (CairoMakie):** a `NaN` or zero-width `colorrange` crashes the PDF renderer; it is built from finite values only and widened when degenerate, with the
  degeneracy check in Float32 (CairoMakie casts down; two close Float64 values can round to the same Float32). Grouped legends with titles cannot be saved to PDF,
  and Makie's nonlinear axis scales make dashed lines ragged, so `optimizer_tests.jl` transforms the data itself (symlog in y, log10 in x) and uses linear axes.
* **Memory:** each worker needs about 2 GB for a DeepCore-sized likelihood; running several multi-worker jobs at once can end in `OutOfMemoryError` inside the
  likelihood evaluation (BAT wraps it in a confusing `MethodError`).
* **Windows:** `time()` has a coarse resolution, so very fast fits can report 0 s.
