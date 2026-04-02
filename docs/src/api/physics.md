# Physics API

## Oscillations

The oscillation module is accessed via `Newtrinos.osc`. Key types include `OscillationConfig`, flavour models (`ThreeFlavour`, `Sterile`, `ADD`), interaction models (`Vacuum`, `SI`, `NSI`), and propagation models (`Basic`, `Decoherent`, `Damping`).

'''@autodocs
Modules = [Newtrinos.osc]
'''

## Eigendecomposition

```@docs
Newtrinos.BargerEigen
```

## Atmospheric Flux

The atmospheric flux module is accessed via `Newtrinos.atm_flux`.

## Cross-Sections

The cross-section module is accessed via `Newtrinos.xsec`.

## Utility Functions

```@docs
Newtrinos.Helpers.bin
Newtrinos.Helpers.rebin
```
