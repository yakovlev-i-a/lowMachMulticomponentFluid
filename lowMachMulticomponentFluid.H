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

Class
    Foam::solvers::lowMachMulticomponentFluid

Description
    Solver module for transient laminar/turbulent flow of a multicomponent
    reacting fluid in the low-Mach approximation.

    Differences from multicomponentFluid:
      - the thermodynamic pressure p is held constant at pRef;
      - only the dynamic pressure pd is solved (the acoustic term psi*ddt(p)
        is removed, i.e. acoustics is filtered out);
      - the density rho = rho(pRef, T, Y) does not depend on pd.

    This is analogous to the "incompressible ideal gas" option in Fluent and
    the p0/pd pressure splitting of lowMachBuoyantPimpleFoam.

    Equations:
        ddt(rho) + div(phiHbyA) - laplacian(rhorAUf, pd) = sources
        U = HbyA + rAU*reconstruct((phi - phiHbyA)/rhorAUf) [+ buoyancy]

SourceFiles
    lowMachMulticomponentFluid.C

See also
    Foam::solvers::multicomponentFluid

\*---------------------------------------------------------------------------*/

#ifndef lowMachMulticomponentFluid_H
#define lowMachMulticomponentFluid_H

#include "multicomponentFluid.H"
#include "pressureReference.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                    Class lowMachMulticomponentFluid Declaration
\*---------------------------------------------------------------------------*/

class lowMachMulticomponentFluid
:
    public multicomponentFluid
{

protected:

    // Dynamic (mechanical) pressure [Pa]; the actually-solved field.
    volScalarField pd_;

    // Thermodynamic reference pressure pRef [Pa] (constant; read from
    // constant/pRef; defaults to 0).
    uniformDimensionedScalarField pRef_;

    // Reference cell/level control for the pd equation.
    Foam::pressureReference pdReference_;


    // Private Member Functions

        //- Construct and solve the low-Mach dynamic pressure equation and
        //  correct the velocity.
        void correctDynamicPressure();


public:

    // Public Data

        //- Reference to the dynamic pressure.
        const volScalarField& pd;


    //- Runtime type information
    TypeName("lowMachMulticomponentFluid");


    // Constructors

        //- Construct from region mesh
        lowMachMulticomponentFluid(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        lowMachMulticomponentFluid(const lowMachMulticomponentFluid&) = delete;


    //- Destructor
    virtual ~lowMachMulticomponentFluid();


    // Member Functions

        //- Called at the start of the PIMPLE loop:
        //  hold the thermodynamic pressure at pRef
        //  (the density is updated inside pressureCorrector).
        virtual void prePredictor();

        //- Predict momentum using the dynamic pressure pd.
        virtual void momentumPredictor();

        //- Solve the low-Mach dynamic pressure equation and correct velocity.
        virtual void pressureCorrector();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const lowMachMulticomponentFluid&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
