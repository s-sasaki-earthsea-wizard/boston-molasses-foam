/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
    Copyright (C) 2026 Syota Sasaki
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    interThermalFoam

Group
    grpMultiphaseSolvers

Description
    Solver for two incompressible, immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    extended with

    - a mixture energy equation in temperature form

        d(rho*Cp*T)/dt + div(rho*Cp*U*T) = div(kappa*grad(T))

      with rho*Cp and kappa volume-averaged by the phase fraction, and

    - an Andrade/Arrhenius temperature-dependent viscosity of phase 1

        mu1(T) = mu0*exp(B*(1/T - 1/T0))
        muEff  = alpha1*mu1(T) + (1 - alpha1)*mu2
        nuEff  = muEff/rho

    Derived from interFoam (OpenFOAM v2312). Laminar flow, static mesh.

    Model switches and coefficients are read from
    constant/transportProperties:
    \verbatim
        water { ...; Cp <J/kg/K>; kappa <W/m/K>; }
        air   { ...; Cp <J/kg/K>; kappa <W/m/K>; }

        energy
        {
            solveEnergy          yes;
        }

        liquidViscosity
        {
            temperatureDependent yes;
            mu0                  5;       // [Pa s] at T0
            T0                   293.15;  // [K]
            B                    1e4;     // [K]
            TMin                 250;     // [K]  (optional)
            TMax                 373.15;  // [K]  (optional)
            muMax                1e4;     // [Pa s] (optional)
        }
    \endverbatim

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "pimpleControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, immiscible fluids using VOF"
        " phase-fraction based interface capturing, with a mixture energy"
        " equation and a temperature-dependent (Andrade) viscosity of"
        " phase 1. Laminar, static mesh."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "initContinuityErrs.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createAlphaFluxes.H"

    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture.correct();

            if (pimple.frozenFlow())
            {
                continue;
            }

            // Energy first, so that the momentum equation sees mu(T^{n+1})
            #include "TEqn.H"
            #include "updateViscosity.H"

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
