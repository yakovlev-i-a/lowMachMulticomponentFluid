/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2026 (low-Mach fork)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

\*---------------------------------------------------------------------------*/

#include "lowMachMulticomponentFluid.H"
#include "addToRunTimeSelectionTable.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "fvcFlux.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"
#include "fvcMeshPhi.H"
#include "fvmDiv.H"
#include "fvmDdt.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(lowMachMulticomponentFluid, 0);
    addToRunTimeSelectionTable(solver, lowMachMulticomponentFluid, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::lowMachMulticomponentFluid::lowMachMulticomponentFluid
(
    fvMesh& mesh
)
:
    multicomponentFluid(mesh),

    pd_
    (
        IOobject
        (
            "pd",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    pRef_
    (
        IOobject
        (
            "pRef",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar(dimPressure, 0)
    ),

    pdReference_(pd_, pimple.dict(), true),

    pd(pd_)
{
    // Hold the thermodynamic pressure at the constant reference value pRef.
    p_ = pRef_;
    p_.correctBoundaryConditions();

    // Initial density from the equation of state at the constant p = pRef.
    rho_ = thermo_.rho();

    // The dynamic pressure pd requires a flux for the pressure-corrector
    // flux reconstruction.
    mesh.schemes().setFluxRequired(pd_.name());

    Info<< "lowMachMulticomponentFluid: pRef = " << pRef_.value() << nl
        << "lowMachMulticomponentFluid: thermodynamic pressure held at pRef, "
        << "solving dynamic pressure pd" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::lowMachMulticomponentFluid::~lowMachMulticomponentFluid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::lowMachMulticomponentFluid::prePredictor()
{
    // Low-Mach approximation: the thermodynamic pressure stays constant.
    // The density is updated inside correctDynamicPressure() from thermo.rho()
    // at p = pRef, so no separate density transport equation is solved here.
    p_ = pRef_;
    p_.correctBoundaryConditions();
}


void Foam::solvers::lowMachMulticomponentFluid::momentumPredictor()
{
    volVectorField& U(U_);

    tUEqn =
    (
        fvm::ddt(rho_, U) + fvm::div(phi_, U)
      + MRF.DDt(rho_, U)
      + momentumTransport->divDevTau(U)
     ==
        fvModels().source(rho_, U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvConstraints().constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        if (buoyancy.valid())
        {
            // Buoyancy force plus the dynamic-pressure gradient
            // (cf. lowMachBuoyantPimpleFoam).
            solve
            (
                UEqn
             ==
                fvc::reconstruct
                (
                    (-buoyancy->ghf*fvc::snGrad(rho_) - fvc::snGrad(pd_))
                   *mesh.magSf()
                )
            );
        }
        else
        {
            solve(UEqn == -fvc::grad(pd_));
        }

        fvConstraints().constrain(U);
        K = 0.5*magSqr(U);
    }
}


void Foam::solvers::lowMachMulticomponentFluid::correctDynamicPressure()
{
    volScalarField& rho(rho_);
    volScalarField& pd(pd_);
    volVectorField& U(U_);
    surfaceScalarField& phi(phi_);

    // Density from the equation of state at the constant p = pRef.
    rho = thermo_.rho();
    rho.relax();

    fvVectorMatrix& UEqn = tUEqn.ref();

    const surfaceScalarField rhof(fvc::interpolate(rho));

    const volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));

    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, pd));

    if (pimple.nCorrPiso() <= 1)
    {
        tUEqn.clear();
    }

    // Standard OpenFOAM mass-flux interpolation:
    // interpolate(rho) * interpolate(HbyA).
    const surfaceScalarField phiHbyAFlux
    (
        "phiHbyAFlux",
        rhof*fvc::flux(HbyA)
    );

    const surfaceScalarField phiHbyADdtCorr
    (
        "phiHbyADdtCorr",
        rhorAUf*fvc::ddtCorr(rho, U, phi, rhoUf)
    );

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        phiHbyAFlux + phiHbyADdtCorr
    );

    MRF.makeRelative(rhof, phiHbyA);

    tmp<surfaceScalarField> tphig;
    if (buoyancy.valid())
    {
        // Buoyancy flux (kept separately for the momentum corrector).
        tphig = - rhorAUf*buoyancy->ghf*fvc::snGrad(rho)*mesh.magSf();
        phiHbyA += tphig();
    }

    constrainPressure(pd, rho, U, phiHbyA, rhorAUf, MRF);
    fvc::makeRelative(phiHbyA, rho, U);

    // Low-Mach dynamic pressure equation (strict Poisson): the acoustic term
    // psi*ddt(p) is removed and the density comes from the EOS at p = pRef.
    fvScalarMatrix pdEqn
    (
        fvc::ddt(rho)
      + fvc::div(phiHbyA)
      - fvm::laplacian(rhorAUf, pd)
     ==
        fvModels().sourceProxy(rho, pd)
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pdEqnN(pdEqn);

        pdEqnN.setReference
        (
            pdReference_.refCell(),
            pdReference_.refValue()
        );

        fvConstraints().constrain(pdEqnN);

        pdEqnN.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pdEqnN.flux();
        }
    }

    // Explicit relaxation of the dynamic pressure for the momentum corrector.
    pd.relax();

    // Solve the continuity equation for the density so that rho remains
    // discretely consistent with the mass flux phi (analogue of the
    // correctDensity() step in the compressible isothermalFluid solver).
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho) + fvc::div(phi)
     ==
        fvModels().source(rho)
    );

    fvConstraints().constrain(rhoEqn);

    rhoEqn.solve();

    fvConstraints().constrain(rho);

    // Velocity correction: reconstruct from the face mass flux for discrete
    // consistency with the mass flux under variable density (v1812
    // lowMachBuoyantPimpleFoam). Using fvc::grad(pd) instead breaks the
    // mass-flux consistency across a sharp density front and smears the
    // combustion products.
    if (buoyancy.valid())
    {
        U = HbyA + rAU*fvc::reconstruct((phi - phiHbyA + tphig())/rhorAUf);
    }
    else
    {
        U = HbyA + rAU*fvc::reconstruct((phi - phiHbyA)/rhorAUf);
    }

    U.correctBoundaryConditions();
    fvConstraints().constrain(U);
    K = 0.5*magSqr(U);
}


void Foam::solvers::lowMachMulticomponentFluid::pressureCorrector()
{
    while (pimple.correct())
    {
        correctDynamicPressure();
    }

    tUEqn.clear();
}


// ************************************************************************* //
