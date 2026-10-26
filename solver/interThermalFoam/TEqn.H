// Mixture energy equation in temperature form:
//
//     ddt(rhoCp, T) + div(rhoCpPhi, T) = laplacian(kappaEff, T)
//
// rhoCp and kappaEff are volume (phase-fraction) averages. The heat-capacity
// flux is rebuilt from the MULES mass flux so that it is consistent with the
// interface advection, including alpha sub-cycling:
//
//     alphaPhi1 = (rhoPhi - rho2*phi)/(rho1 - rho2)
//     rhoCpPhi  = alphaPhi1*(rhoCp1 - rhoCp2) + rhoCp2*phi
//
// The residual of the rhoCp continuity equation is removed implicitly, which
// turns the operator into rhoCp*DT/Dt and keeps T bounded.

if (solveEnergy)
{
    // Conductivity must stay positive: use the bounded phase fraction
    const volScalarField limitedAlpha1
    (
        "limitedAlpha1",
        min(max(alpha1, scalar(0)), scalar(1))
    );

    rhoCp == alpha1*rhoCp1 + alpha2*rhoCp2;

    kappaEff == limitedAlpha1*kappa1 + (scalar(1) - limitedAlpha1)*kappa2;

    rhoCpPhi =
        (rhoPhi - rho2*phi)*((rhoCp1 - rhoCp2)/(rho1 - rho2))
      + rhoCp2*phi;

    fvScalarMatrix TEqn
    (
        fvm::ddt(rhoCp, T)
      + fvm::div(rhoCpPhi, T)
      - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoCpPhi), T)
      - fvm::laplacian(kappaEff, T)
     ==
        fvOptions(rhoCp, T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);

    Info<< "T min/max = " << gMin(T.primitiveField())
        << ", " << gMax(T.primitiveField()) << endl;
}
