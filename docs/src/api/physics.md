# Physics API

## Oscillations

The oscillation module is accessed via `Newtrinos.osc`. Use [`Newtrinos.osc.configure`](@ref) to
create an [`Newtrinos.osc.Osc`](@ref) physics module that can be passed to any experiment.

### Configuration

```@docs
Newtrinos.osc.OscillationConfig
Newtrinos.osc.Osc
Newtrinos.osc.configure
```

### Flavour Models

```@docs
Newtrinos.osc.FlavourModel
Newtrinos.osc.ThreeFlavour
Newtrinos.osc.ThreeFlavourXYCP
Newtrinos.osc.Sterile
Newtrinos.osc.ADD
```

### Interaction Models

```@docs
Newtrinos.osc.InteractionModel
Newtrinos.osc.Vacuum
Newtrinos.osc.SI
Newtrinos.osc.NSI
```

### Propagation Models

```@docs
Newtrinos.osc.PropagationModel
Newtrinos.osc.Basic
Newtrinos.osc.Decoherent
Newtrinos.osc.Damping
```

### State Selection

```@docs
Newtrinos.osc.StateSelector
Newtrinos.osc.All
Newtrinos.osc.Cut
```

### Eigendecomposition

```@docs
Newtrinos.osc.EigenMethod
Newtrinos.osc.DefaultEigen
Newtrinos.osc.decompose
Newtrinos.BargerEigen
```

### Data Types

```@docs
Newtrinos.osc.ftype
Newtrinos.osc.Layer
Newtrinos.osc.Path
```

### Core Functions

```@docs
Newtrinos.osc.get_osc_prob
Newtrinos.osc.get_matrices
Newtrinos.osc.get_PMNS
Newtrinos.osc.get_abs_masses
```

### Parameters and Priors

```@docs
Newtrinos.osc.get_params
Newtrinos.osc.get_priors
```

### Internal Functions

These functions are not part of the public API but are documented for developers working
on the oscillation engine.

```@docs
Newtrinos.osc.osc_kernel
Newtrinos.osc.compute_matter_matrices
Newtrinos.osc.osc_reduce
Newtrinos.osc.matter_osc_per_e
Newtrinos.osc.select
Newtrinos.osc.propagate
```

## Atmospheric Flux

The atmospheric flux module is accessed via `Newtrinos.atm_flux`.

## Cross-Sections

The cross-section module is accessed via `Newtrinos.xsec`. Use [`Newtrinos.xsec.configure`](@ref)
to create a [`Newtrinos.xsec.Xsec`](@ref) physics module.

### Models

```@docs
Newtrinos.xsec.XsecModel
Newtrinos.xsec.SimpleScaling
Newtrinos.xsec.Differential_H2O
Newtrinos.xsec.Xsec
```

### Configuration and Functions

```@docs
Newtrinos.xsec.configure
Newtrinos.xsec.get_params
Newtrinos.xsec.get_priors
Newtrinos.xsec.get_scale
```

## Utility Functions

```@docs
Newtrinos.Helpers.bin
Newtrinos.Helpers.rebin
```
