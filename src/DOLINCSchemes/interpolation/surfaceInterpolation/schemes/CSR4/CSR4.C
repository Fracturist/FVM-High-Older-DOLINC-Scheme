/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

Description
    4th-order central-stencil reconstruction scheme class

Author
    Hao Guo

\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "CSR4.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::surfaceScalarField>
Foam::CSR4<Type>::weights
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    const surfaceScalarField& CDWeights = mesh.weights();

    tmp<surfaceScalarField> tsWeights
    (
        new surfaceScalarField
        (
            IOobject
            (
                "CSR4::schemeWeights(" + vf.name() + ')',
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false  // Do not register
            ),
            mesh,
            dimless
        )
    );

    surfaceScalarField& sWeights = tsWeights.ref();
    sWeights.setOriented();

    forAll(sWeights, faceI)
    {
        scalar g = CDWeights[faceI]/(1.0 - CDWeights[faceI]);

        sWeights[faceI] = schemeWeights(g);
    }

    surfaceScalarField::Boundary& bSWeights = sWeights.boundaryFieldRef();

    forAll(bSWeights, patchI)
    {
        fvsPatchField<scalar>& pSWeights = bSWeights[patchI];

        if (pSWeights.coupled())
        {
            const scalarField& pCDWeights = CDWeights.boundaryField()[patchI];

            forAll(pSWeights, faceI)
            {
                scalar g = pCDWeights[faceI]/(1.0 - pCDWeights[faceI]);

                pSWeights[faceI] = schemeWeights(g);
            }
        }
        else
        {
            pSWeights = 1.0;
        }
    }

    return tsWeights;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::CSR4<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    const surfaceScalarField& CDWeights = mesh.weights();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsfCorr
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "CSR4::correction(" + vf.name() + ')',
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensioned<Type>(vf.name(), vf.dimensions(), Zero)
        )
    );

    GeometricField<Type, fvsPatchField, surfaceMesh>& sfCorr = tsfCorr.ref();

    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();

    tmp<fv::gradScheme<scalar>> gradScheme_
    (
        fv::gradScheme<scalar>::New
        (
            mesh,
            mesh.gradScheme(gradSchemeName_)
        )
    );

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        tmp<volVectorField> tgradVf =
            gradScheme_().grad(vf.component(cmpt), gradSchemeName_);

        const volVectorField& gradVf = tgradVf();

        forAll(sfCorr, faceI)
        {
            const label cellO = owner[faceI];
            const label cellN = neighbour[faceI];

            scalar g = CDWeights[faceI]/(1.0 - CDWeights[faceI]);

            setComponent(sfCorr[faceI], cmpt) =
                schemeCorr
                (
                    g,
                    gradVf[cellO],
                    gradVf[cellN],
                    C[cellN] - C[cellO]
                );
        }

        typename GeometricField<Type, fvsPatchField, surfaceMesh>::
            Boundary& bSfCorr = sfCorr.boundaryFieldRef();

        forAll(bSfCorr, patchI)
        {
            fvsPatchField<Type>& pSfCorr = bSfCorr[patchI];

            if (pSfCorr.coupled())
            {
                const scalarField& pCDWeights = CDWeights.boundaryField()[patchI];

                const labelUList& pOwner = mesh.boundary()[patchI].faceCells();

                const vectorField pGradVfNei
                (
                    gradVf.boundaryField()[patchI].patchNeighbourField()
                );

                // Build the d-vectors
                const vectorField pd
                (
                    vf.boundaryField()[patchI].patch().delta()
                );

                forAll(pOwner, faceI)
                {
                    // label own = pOwner[faceI];

                    scalar g = pCDWeights[faceI]/(1.0 - pCDWeights[faceI]);

                    setComponent(pSfCorr[faceI], cmpt) =
                        schemeCorr
                        (
                            g,
                            gradVf[pOwner[faceI]],
                            pGradVfNei[faceI],
                            pd[faceI]
                        );
                }
            }
            else
            {
                forAll(pSfCorr, faceI)
                {
                    setComponent(pSfCorr[faceI], cmpt) = 0.0;
                }
            }
        }
    }

    return tsfCorr;
}


namespace Foam
{
    makeSurfaceInterpolationScheme(CSR4)
}

// ************************************************************************* //
