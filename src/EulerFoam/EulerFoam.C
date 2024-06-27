/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Application
    EulerFoam

Group
    grpCompressibleSolvers

Description
    Density-based compressible flow solver based on central-upwind
    schemes of Kurganov and Tadmor.

Author
    Hao Guo

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
#include "directionInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Density-based compressible flow solver based on central-upwind"
        " schemes of Kurganov and Tadmor."
    );

    argList::addBoolOption
    (
        "debugFields",
        "Write extra fields for debugging"
    );

    #define NO_CONTROL
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFieldRefs.H"

    #include "createTimeControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "readAdvanceScheme.H"
    #include "readFluxScheme.H"

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        surfaceScalarField fluxRho
        (
            IOobject
            (
                "flux(rho)",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar(dimDensity*dimVelocity*dimArea, Zero)
        );

        surfaceVectorField fluxRhoU
        (
            IOobject
            (
                "flux(rhoU)",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedVector(dimPressure*dimArea, Zero)
        );

        surfaceScalarField fluxRhoE
        (
            IOobject
            (
                "flux(rhoE)",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar(dimPressure*dimVelocity*dimArea, Zero)
        );

        #include "defineAndReconstructRiemannVars.H"
        #include "EulerCourantNo.H"

        // For adjustable time step
        #include "readTimeControls.H"  // update adjustable time step settings
        #include "setDeltaT.H"  // need CoNum

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (advanceScheme == "ForwardEuler")
        {
            #include "E1.H"
        }
        else if (advanceScheme == "TVDRungeKutta3")
        {
            #include "TVDRK3.H"
        }

        if (debugF)
        {
            e.write();

            psi.write();
            rhoU.write();
            rhoE.write();

            fluxRho.write();
            fluxRhoU.write();
            fluxRhoE.write();
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
