/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Author
    Hao Guo

\*---------------------------------------------------------------------------*/

#include "firstOrderDolincScheme.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "coupledFvPatchField.H"


// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class Type>
void Foam::firstOrderDolincScheme<Type>::calcDiff
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    faceOwnDiff1_.reset
    (
        autoPtr<GeometricField<Type, fvsPatchField, surfaceMesh>>::New
        (
            IOobject
            (
                "firstOrderDolincScheme::owner1stOrderDiff(" + vf.name() + ")",
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

    faceNeiDiff1_.reset
    (
        autoPtr<GeometricField<Type, fvsPatchField, surfaceMesh>>::New
        (
            IOobject
            (
                "firstOrderDolincScheme::neighbour1stOrderDiff(" + vf.name() + ")",
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

    GeometricField<Type, fvsPatchField, surfaceMesh>& faceOwnDiff1 = faceOwnDiff1_();
    GeometricField<Type, fvsPatchField, surfaceMesh>& faceNeiDiff1 = faceNeiDiff1_();

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

        forAll(faceOwnDiff1, faceI)
        {
            label fOwner = owner[faceI];
            label fNeighbour= neighbour[faceI];

            const vector twoDelta = 2*(C[fNeighbour] - C[fOwner]);

            setComponent(faceOwnDiff1[faceI], cmpt) =
                gradVf[fOwner] & twoDelta;
            setComponent(faceNeiDiff1[faceI], cmpt) =
                gradVf[fNeighbour] & twoDelta;
        }

        typename GeometricField<Type, fvsPatchField, surfaceMesh>::
            Boundary& bFaceOwnDiff1 = faceOwnDiff1.boundaryFieldRef();
        typename GeometricField<Type, fvsPatchField, surfaceMesh>::
            Boundary& bFaceNeiDiff1 = faceNeiDiff1.boundaryFieldRef();

        forAll(bFaceOwnDiff1, patchI)
        {
            fvsPatchField<Type>& pFaceOwnDiff1 = bFaceOwnDiff1[patchI];

            if (pFaceOwnDiff1.coupled())
            {
                fvsPatchField<Type>& pFaceNeiDiff1 = bFaceNeiDiff1[patchI];

                const labelUList& pOwner = mesh.boundary()[patchI].faceCells();

                const vectorField pGradVfNei
                (
                    gradVf.boundaryField()[patchI].patchNeighbourField()
                );

                // Build the patch 2d vectors
                const vectorField p2d
                (
                    2*vf.boundaryField()[patchI].patch().delta()
                );

                forAll(pOwner, faceI)
                {
                    label owner = pOwner[faceI];

                    setComponent(pFaceOwnDiff1[faceI], cmpt) =
                        gradVf[owner] & p2d[faceI];
                    setComponent(pFaceNeiDiff1[faceI], cmpt) =
                        pGradVfNei[faceI] & p2d[faceI];
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::firstOrderDolincScheme<Type>>
Foam::firstOrderDolincScheme<Type>::New
(
    const fvMesh& mesh,
    Istream& schemeData
)
{
    if (surfaceInterpolation::debug)
    {
        InfoInFunction
            << "Constructing firstOrderDolincScheme<Type>" << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorInFunction(schemeData)
            << "Discretisation scheme not specified"
            << endl << endl
            << "Valid schemes are :" << endl
            << MeshConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    auto cstrIter = MeshConstructorTablePtr_->cfind(schemeName);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            schemeData,
            "discretisation",
            schemeName,
            *MeshConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return cstrIter()(mesh, schemeData);
}


template<class Type>
Foam::tmp<Foam::firstOrderDolincScheme<Type>>
Foam::firstOrderDolincScheme<Type>::New
(
    const fvMesh& mesh,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
)
{
    if (surfaceInterpolation::debug)
    {
        InfoInFunction
            << "Constructing firstOrderDolincScheme<Type>"
            << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorInFunction(schemeData)
            << "Discretisation scheme not specified"
            << endl << endl
            << "Valid schemes are :" << endl
            << MeshConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    auto cstrIter = MeshFluxConstructorTablePtr_->cfind(schemeName);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            schemeData,
            "discretisation",
            schemeName,
            *MeshFluxConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return cstrIter()(mesh, faceFlux, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::firstOrderDolincScheme<Type>::
~firstOrderDolincScheme()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// template<class Type>
// Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
// Foam::firstOrderDolincScheme<Type>::flux
// (
//     const GeometricField<Type, fvPatchField, volMesh>& phi
// ) const
// {
//     return faceFlux_*this->interpolate(phi);
// }


// ************************************************************************* //
