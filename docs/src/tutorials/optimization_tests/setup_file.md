# Creating a setup file

The analysis scripts (`optimizer_tests.jl`, `roundtrips.jl`, `nsi/nsi_scan.jl`) contain no physics. They are given a **setup file** with `--setup`, which
describes the *problem*: which experiments and physics model, which parameters and priors, and which questions (hypotheses) to ask. How the problem is
solved (optimizers, budgets, seeds) stays on the command line. To analyse something else (another dataset, another model) you only write a new setup file.

A setup file is an ordinary Julia file that defines one function:

```julia
setup(; ftype = Float64) -> (; name, experiments, params, priors, hypotheses)
```

| field | content |
|---|---|
| `name` | short label; used for the default results folder `results/<name>/...` and stored in every saved result |
| `experiments` | NamedTuple of configured experiments, all built from the same `physics` |
| `params` | start values of all parameters, `Newtrinos.get_params(experiments)` |
| `priors` | priors of all parameters (before any hypothesis is applied), `Newtrinos.get_priors(experiments)` |
| `hypotheses` | ordered dictionary `name => (; prior, null)`, see step 4 |

The scripts `include` the file and call `setup(; ftype = Float64)`, so the file may load packages and define helper functions at the top level.

## Step 1: physics and experiments

Configure the physics and the experiments exactly as in [Configure an experiment](../configure_experiment.md), and keep them in the NamedTuple `experiments`:

```julia
using Distributions, DataStructures, Accessors, Newtrinos

osc = Newtrinos.osc.configure()                                   # three flavours, vacuum
experiments = (; dayabay = Newtrinos.dayabay.configure((; osc)))
```

For a joint analysis list several experiments in the tuple; the likelihood is the product (see the *Joint Analysis* example).

## Step 2: parameters and priors

```julia
p      = Newtrinos.get_params(experiments)
priors = Newtrinos.get_priors(experiments)
```

`priors` holds one distribution per parameter. Two operations adapt it, both applied to `priors` **before** the hypotheses are built (every hypothesis is
derived from `priors`):

* **Fix a parameter** (it becomes a constant in every fit and every truth): `priors = Newtrinos.condition(priors, Dict(:θ₁₂ => p.θ₁₂), params)`.
* **Change a prior** (the parameter stays free): `@reset priors.θ₂₃ = Uniform(ftype(0.6), ftype(1.0))`.

Fixing is how you keep parameters that the experiment is not sensitive to (or that are not part of the question) out of the optimization.

## Step 3: what is a hypothesis

A hypothesis is one question of the form *"is effect X present, and can the fit recover it?"*. It has two parts, and both are a plain NamedTuple field:

* **`prior`**: the prior in which the parameters of X (the **tested** parameters) are **free**, together with the nuisance parameters you want fitted, while
  everything else is **fixed** (`Newtrinos.condition`).
* **`null`**: a NamedTuple with the **tested** parameters and the values they take when the effect X is absent, e.g. `(ε_eμ_abs = 0.0, δ_eμ = 0.0)` for
  "no non-standard interaction", or the nominal value for a parameter that is always present, e.g. `(θ₁₃ = p.θ₁₃,)`.

Everything else is derived by the scripts:

| derived | rule |
|---|---|
| tested parameters | `keys(null)`; their order is the order of the scan axes (2D: first name = x axis) |
| scan type | 1 tested parameter → 1D grid, 2 → 2D grid, more → Sobol sample |
| null model | `condition(prior, null, params)`: the same prior with X switched off; used by the `smseed` stage |
| nuisance parameters | all other free parameters of `prior`; they are fitted at every point and redrawn from their priors by `--fluctuate-nuisances` |

## Step 4: build the hypotheses

Write a small helper that builds a hypothesis from the list of tested parameters and the list of parameters that stay free:

```julia
function hyp(tested; free)
    pinned = Dict{Symbol, Any}(k => p[k] for k in (:θ₁₃, :θ₂₃, :Δm²₃₁) if !(k in free))   # fix everything that is not free
    prior  = isempty(pinned) ? deepcopy(priors) : Newtrinos.condition(priors, pinned, p)
    (; prior, null = NamedTuple{Tuple(tested)}(Tuple(p[k] for k in tested)))
end

hypotheses = OrderedDict(                                              # ordered: the scripts list and run hypotheses in this order
    "theta13"      => hyp((:θ₁₃,);             free = (:θ₁₃, :Δm²₃₁)),   # 1D scan; Δm²₃₁ is a nuisance parameter
    "theta13_dm31" => hyp((:θ₁₃, :Δm²₃₁);      free = (:θ₁₃, :Δm²₃₁)),   # 2D scan
    "three"        => hyp((:θ₁₃, :Δm²₃₁, :θ₂₃); free = (:θ₁₃, :θ₂₃, :Δm²₃₁)),   # Sobol sample
)
```

The last line of the file's function returns the NamedTuple:

```julia
(; name = "example_dayabay", experiments, params = p, priors, hypotheses)
```

## The complete example

```julia
# my_setup.jl: Daya Bay with three-flavour vacuum oscillations
using Distributions, DataStructures, Accessors, Newtrinos

function setup(; ftype = Float64)
    osc = Newtrinos.osc.configure()
    experiments = (; dayabay = Newtrinos.dayabay.configure((; osc)))

    p      = Newtrinos.get_params(experiments)
    priors = Newtrinos.get_priors(experiments)

    # Daya Bay is insensitive to these: fixed in every hypothesis (δCP stays free: a nuisance parameter)
    priors = Newtrinos.condition(priors, Dict(:θ₁₂ => p.θ₁₂, :Δm²₂₁ => p.Δm²₂₁), p)

    function hyp(tested; free)
        pinned = Dict{Symbol, Any}(k => p[k] for k in (:θ₁₃, :θ₂₃, :Δm²₃₁) if !(k in free))
        prior  = isempty(pinned) ? deepcopy(priors) : Newtrinos.condition(priors, pinned, p)
        (; prior, null = NamedTuple{Tuple(tested)}(Tuple(p[k] for k in tested)))
    end

    hypotheses = OrderedDict(
        "theta13"      => hyp((:θ₁₃,);              free = (:θ₁₃, :Δm²₃₁)),
        "theta13_dm31" => hyp((:θ₁₃, :Δm²₃₁);       free = (:θ₁₃, :Δm²₃₁)),
        "three"        => hyp((:θ₁₃, :Δm²₃₁, :θ₂₃); free = (:θ₁₃, :θ₂₃, :Δm²₃₁)),
    )
    (; name = "my_setup", experiments, params = p, priors, hypotheses)
end
```

Use it with any of the scripts:

```bash
julia src/analysis/optimizer_tests.jl --setup my_setup.jl --hypothesis theta13 --strategies lbfgsb "mh-de lbfgsb" --workers 3
julia src/analysis/roundtrips.jl      --setup my_setup.jl --hypothesis theta13_dm31 --strategy random lbfgsb --n-points 16
julia nsi/nsi_scan.jl                 --setup my_setup.jl --name run1 --task profile --hypothesis theta13 --objective posterior --strategy random lbfgsb
```

A larger real setup with many hypotheses, detector systematics with their paper priors and a hypothesis helper that pins all other parameters to their null value
is `nsi/setups/nsi_deepcore_3y.jl`.

## Rules and common mistakes

The scripts check the setup (`finalize_setup`) and stop with a message if a rule is broken. These are the ones you will meet:

* **A tested parameter must be free in the hypothesis prior.** If the base `priors` fixes a parameter (step 2), no hypothesis can test it; the scripts stop with
  *"tested parameter θ₂₃ is fixed in its prior, it must be free"*. Fix parameters in the base `priors` only if *no* hypothesis tests them; pin the others
  per hypothesis (as `hyp` does above).
* **Keep at least one nuisance parameter free** in every hypothesis you scan with a *profile*: at each grid point the tested parameters are pinned, and if nothing
  else is free there is nothing left to fit. The point is then just evaluated, without an optimization, which is usually not what you want. (The 2D and the
  3-parameter hypotheses above stay fittable only because δCP is left free.)
* **The null value must be a sensible point of the prior**, typically the value of the parameter without the effect, or its nominal value. It is used for the null model of
  `smseed`, and is the "effect absent" reference.
* **Do the fixing and prior changes on `priors` first**, then build the hypotheses: a hypothesis is a copy of `priors` with extra parameters pinned.
* **Use `ftype` for literals** in priors (`Uniform(ftype(0), ftype(2π))`) if the file is meant to run in other precisions.
* **Hypothesis names** are what `--hypothesis` takes; non-ASCII names such as `"ττ_μμ"` are fine. Without `--hypothesis`, `optimizer_tests.jl` and
  `roundtrips.jl` use the first hypothesis of the dictionary.

## Checking a setup quickly

Before a long run, load the file and look at what the scripts will derive from it:

```julia
using Distributions, DensityInterface, BAT, DataStructures, MeasureBase, ADTypes, ValueShapes, Accessors, Newtrinos
import Optim, ForwardDiff                              # the packages the scripts load before the helper file
include("src/analysis/optimizer_tests_common.jl")     # the helpers of the scripts
include("my_setup.jl")
S = finalize_setup(setup())                           # checks the fields and every hypothesis
for name in keys(S.hypotheses)
    h = hypothesis_info(S, name)
    println(name, ": tested = ", h.tested, ", nuisance = ", h.fluctuated)
end
likelihood = Newtrinos.generate_likelihood(S.experiments)
DensityInterface.logdensityof(likelihood, S.params)   # the likelihood at the start values must be finite
```
