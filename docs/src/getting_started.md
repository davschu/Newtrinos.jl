# Getting Started

## Installation

Newtrinos.jl is not yet registered in the Julia General registry. Install it directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/philippeller/Newtrinos.jl.git")
```

Or for development:

```bash
git clone https://github.com/philippeller/Newtrinos.jl.git
cd Newtrinos.jl
julia --project -e 'using Pkg; Pkg.instantiate()'
```

## First Steps

Load the package:

```julia
using Newtrinos
using DensityInterface
```

Configure an experiment with default physics:

```julia
exp = Newtrinos.dayabay.configure()
```

Each experiment returns a struct containing physics models, parameters, priors, data assets, and a forward model. Extract parameters and evaluate the likelihood:

```julia
experiments = (; dayabay = exp)
params = Newtrinos.get_params(experiments)
likelihood = Newtrinos.generate_likelihood(experiments)

logdensityof(likelihood, params)
```

## Running Tests

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

## Next Steps

- [Single Experiment](tutorials/single_experiment.md) — configure and explore one experiment in depth
- [Joint Analysis](tutorials/joint_analysis.md) — combine experiments and run scans
- [Architecture](manual/architecture.md) — understand the three-layer design
