using Newtrinos

function configure_physics(ordering::String)
    osc_cfg = Newtrinos.osc.OscillationConfig(
        flavour=Newtrinos.osc.ThreeFlavour(ordering=Symbol(ordering)),
        propagation=Newtrinos.osc.Basic(),
        states=Newtrinos.osc.All(),
        interaction=Newtrinos.osc.SI()
    )
    osc = Newtrinos.osc.configure(osc_cfg)

    atm_flux = Newtrinos.atm_flux.configure(
        Newtrinos.atm_flux.AtmFluxConfig(nominal_model=Newtrinos.atm_flux.HKKM("kam-ally-20-01-mtn-solmin.d"))
    )
    earth_layers = Newtrinos.earth_layers.configure()
    xsec = Newtrinos.xsec.configure(Newtrinos.xsec.Differential_H2O())

    (; osc, atm_flux, earth_layers, xsec)
end

function configure_experiments(experiment_list, physics)
    pairs = (Symbol(lowercase(exp)) => getproperty(getproperty(Newtrinos, Symbol(lowercase(exp))), :configure)(physics) for exp in experiment_list)
    return (; pairs...)
end
