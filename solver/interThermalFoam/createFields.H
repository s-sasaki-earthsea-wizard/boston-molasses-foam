#include "createRDeltaT.H"

if (LTS)
{
    FatalErrorInFunction
        << "Local time stepping (localEuler) is not supported by "
        << "interThermalFoam" << nl
        << exit(FatalError);
}

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

#include "createFvOptions.H"


// * * * * * * * * * * * * * * * Energy equation * * * * * * * * * * * * * //

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

const dictionary& phase1Dict = mixture.subDict(mixture.phase1Name());
const dictionary& phase2Dict = mixture.subDict(mixture.phase2Name());

const dimensionSet dimThermalConductivity(dimPower/dimLength/dimTemperature);

// Specific heat capacities [J/kg/K] and thermal conductivities [W/m/K]
const dimensionedScalar Cp1("Cp", dimSpecificHeatCapacity, phase1Dict);
const dimensionedScalar Cp2("Cp", dimSpecificHeatCapacity, phase2Dict);
const dimensionedScalar kappa1("kappa", dimThermalConductivity, phase1Dict);
const dimensionedScalar kappa2("kappa", dimThermalConductivity, phase2Dict);

// Volumetric heat capacities [J/m3/K]
const dimensionedScalar rhoCp1("rhoCp1", rho1*Cp1);
const dimensionedScalar rhoCp2("rhoCp2", rho2*Cp2);

const bool solveEnergy
(
    mixture.subDict("energy").get<bool>("solveEnergy")
);

// Mixture volumetric heat capacity; the old-time level is needed by
// fvm::ddt(rhoCp, T)
volScalarField rhoCp
(
    IOobject
    (
        "rhoCp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1*rhoCp1 + alpha2*rhoCp2
);
rhoCp.oldTime();

// Heat-capacity flux [W/K], rebuilt from rhoPhi in TEqn.H
surfaceScalarField rhoCpPhi
(
    IOobject
    (
        "rhoCpPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rhoCp)*phi
);

// Mixture thermal conductivity
volScalarField kappaEff
(
    IOobject
    (
        "kappaEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimThermalConductivity, Zero)
);


// * * * * * * * * * * * * * * * * Viscosity * * * * * * * * * * * * * * * //

const dictionary& viscosityDict = mixture.subDict("liquidViscosity");

const bool temperatureDependentViscosity
(
    viscosityDict.get<bool>("temperatureDependent")
);

// mu1(T) = mu0*exp(B*(1/T - 1/T0))
const dimensionedScalar mu0("mu0", dimDynamicViscosity, viscosityDict);
const dimensionedScalar T0("T0", dimTemperature, viscosityDict);
const dimensionedScalar BVisc("B", dimTemperature, viscosityDict);

// Safety bounds for the exponential law
const dimensionedScalar TMinVisc
(
    "TMin",
    dimTemperature,
    viscosityDict.getOrDefault<scalar>("TMin", 200)
);
const dimensionedScalar TMaxVisc
(
    "TMax",
    dimTemperature,
    viscosityDict.getOrDefault<scalar>("TMax", 500)
);
const dimensionedScalar muMaxVisc
(
    "muMax",
    dimDynamicViscosity,
    viscosityDict.getOrDefault<scalar>("muMax", GREAT)
);

if (mu0.value() <= 0 || T0.value() <= 0 || TMinVisc.value() <= 0)
{
    FatalIOErrorInFunction(viscosityDict)
        << "mu0, T0 and TMin must be positive" << nl
        << exit(FatalIOError);
}

if (temperatureDependentViscosity && !solveEnergy)
{
    WarningInFunction
        << "liquidViscosity/temperatureDependent is on but "
        << "energy/solveEnergy is off: mu(T) is evaluated on the frozen "
        << "initial temperature field" << endl;
}

// The momentum equation uses mu0 / mu(T); the 'nu' entry of phase 1 is only
// read by the two-phase mixture class. Flag inconsistent input.
if (phase1Dict.found("nu"))
{
    const dimensionedScalar nu1("nu", dimViscosity, phase1Dict);
    const scalar muFromNu = (rho1*nu1).value();

    if (mag(muFromNu - mu0.value()) > 1e-3*mu0.value())
    {
        WarningInFunction
            << "rho*nu of phase " << mixture.phase1Name()
            << " (" << muFromNu << " Pa s) differs from liquidViscosity/mu0 ("
            << mu0.value() << " Pa s). The momentum equation uses mu0."
            << endl;
    }
}

Info<< nl << "Thermal/viscosity model:" << nl
    << "    energy equation      : " << (solveEnergy ? "on" : "off") << nl
    << "    liquid viscosity     : "
    << (temperatureDependentViscosity ? "Andrade mu(T)" : "constant") << nl
    << "    mu0, T0, B           : " << mu0.value() << " Pa s, "
    << T0.value() << " K, " << BVisc.value() << " K" << nl
    << "    rho*Cp (" << mixture.phase1Name() << ", "
    << mixture.phase2Name() << ") : " << rhoCp1.value() << ", "
    << rhoCp2.value() << " J/m3/K" << nl
    << "    kappa  (" << mixture.phase1Name() << ", "
    << mixture.phase2Name() << ") : " << kappa1.value() << ", "
    << kappa2.value() << " W/m/K" << nl << endl;

// Phase-1 (liquid) dynamic viscosity
volScalarField muLiquid
(
    IOobject
    (
        "muLiquid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    mu0
);

// Mixture dynamic viscosity used by the momentum equation
volScalarField muEff
(
    IOobject
    (
        "muEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    mu0
);

// Mixture kinematic viscosity (diagnostic)
volScalarField nuEff
(
    IOobject
    (
        "nuEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimViscosity, Zero)
);

#include "updateViscosity.H"
