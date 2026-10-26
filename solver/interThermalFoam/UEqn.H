    // Laminar momentum equation with the variable mixture viscosity muEff
    // (same stress form as the Stokes laminar model's divDevRhoReff).
    // Foam::T() is qualified because the temperature field is also named T.
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      - fvm::laplacian(muEff, U)
      - fvc::div
        (
            muEff*dev2(Foam::T(fvc::grad(U))),
            "div(muEff*dev2(T(grad(U))))"
        )
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
