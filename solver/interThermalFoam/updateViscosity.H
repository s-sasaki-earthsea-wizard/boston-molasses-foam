// Mixture viscosity update:
//
//     mu1(T) = mu0*exp(B*(1/T - 1/T0))            (Andrade/Arrhenius)
//     muEff  = alpha1*mu1(T) + (1 - alpha1)*mu2
//     nuEff  = muEff/rho
//
// As in interFoam's incompressibleTwoPhaseMixture, the phase fraction is
// bounded to [0, 1] for the mixture properties. T is clipped to
// [TMin, TMax] before evaluating the exponential and mu1 is capped at muMax.

{
    const volScalarField limitedAlpha1
    (
        "limitedAlpha1",
        min(max(alpha1, scalar(0)), scalar(1))
    );

    if (temperatureDependentViscosity)
    {
        const volScalarField Tvisc
        (
            "Tvisc",
            min(max(T, TMinVisc), TMaxVisc)
        );

        muLiquid == min(mu0*exp(BVisc*(1.0/Tvisc - 1.0/T0)), muMaxVisc);
    }
    else
    {
        muLiquid == mu0;
    }

    // Phase 2 (air) keeps the viscosity model of the two-phase mixture
    const volScalarField muAir("muAir", rho2*mixture.nuModel2().nu());

    muEff == limitedAlpha1*muLiquid + (scalar(1) - limitedAlpha1)*muAir;

    nuEff ==
        muEff/(limitedAlpha1*rho1 + (scalar(1) - limitedAlpha1)*rho2);
}
