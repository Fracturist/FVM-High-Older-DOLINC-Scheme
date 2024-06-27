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
    2nd order ENO scheme

Author
    Hao Guo

\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "UENO2.H"

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class Type>
void Foam::UENO2<Type>::selectStencil
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    this->checkDiff(vf);

    const fvMesh& mesh = this->mesh();

    const surfaceScalarField& faceFlux = this->faceFlux_;

    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();

    const GeometricField<Type, fvsPatchField, surfaceMesh>& faceOwnDiff1 = this->faceOwnDiff1();
    const GeometricField<Type, fvsPatchField, surfaceMesh>& faceNeiDiff1 = this->faceNeiDiff1();

    stencilNo_.resize(faceFlux.size(), 0);

    forAll(stencilNo_, faceI)
    {
        scalar sm2L1, sm2L0;

        if (faceFlux[faceI] >= 0)
        {
            sm2L1 =
                this->smoothness2L1
                (
                    1.0,
                    vf[owner[faceI]],
                    vf[neighbour[faceI]],
                    faceOwnDiff1[faceI]
                );

            sm2L0 =
                this->smoothness2L0
                (
                    vf[owner[faceI]],
                    vf[neighbour[faceI]]
                );
        }
        else
        {
            sm2L1 =
                this->smoothness2L1
                (
                    1.0,
                    vf[neighbour[faceI]],
                    vf[owner[faceI]],
                    -faceNeiDiff1[faceI]
                );

            sm2L0 =
                this->smoothness2L0
                (
                    vf[neighbour[faceI]],
                    vf[owner[faceI]]
                );
        }

        // stencilNo_[faceI] <= 0
        if ( b_*sm2L1 < sm2L0 )
        {
            stencilNo_[faceI] += 1;
        }
    }

    const typename surfaceScalarField::
        Boundary& bFaceFlux = faceFlux.boundaryField();

    bStencilNo_.resize(bFaceFlux.size(), labelField());

    forAll(bFaceFlux, patchI)
    {
        const fvsPatchField<scalar>& pFaceFlux = bFaceFlux[patchI];

        if (pFaceFlux.coupled())
        {
            labelField& pStencilNo = bStencilNo_[patchI];

            pStencilNo.resize(pFaceFlux.size(), 0);

            const labelUList& pOwner = mesh.boundary()[patchI].faceCells();

            const Field<Type> pVfNei
            (
                vf.boundaryField()[patchI].patchNeighbourField()
            );

            const Field<Type>& pFaceOwnDiff1 = faceOwnDiff1.boundaryField()[patchI];
            const Field<Type>& pFaceNeiDiff1 = faceNeiDiff1.boundaryField()[patchI];

            forAll(pStencilNo, faceI)
            {
                scalar sm2L1, sm2L0;

                if (pFaceFlux[faceI] >= 0)
                {

                    sm2L1 =
                        this->smoothness2L1
                        (
                            1.0,
                            vf[pOwner[faceI]],
                            pVfNei[faceI],
                            pFaceOwnDiff1[faceI]
                        );

                    sm2L0 =
                        this->smoothness2L0
                        (
                            vf[pOwner[faceI]],
                            pVfNei[faceI]
                        );
                }
                else
                {
                    sm2L1 =
                        this->smoothness2L1
                        (
                            1.0,
                            pVfNei[faceI],
                            vf[pOwner[faceI]],
                            -pFaceNeiDiff1[faceI]
                        );

                    sm2L0 =
                        this->smoothness2L0
                        (
                            pVfNei[faceI],
                            vf[pOwner[faceI]]
                        );
                }

                if ( b_*sm2L1 < sm2L0 )
                {
                    pStencilNo[faceI] += 1;
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::surfaceScalarField>
Foam::UENO2<Type>::weights
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
                "UENO2::schemeWeights(" + vf.name() + ')',
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false // Do not register
            ),
            mesh,
            dimless
        )
    );

    surfaceScalarField& sWeights = tsWeights.ref();
    sWeights.setOriented();

    const surfaceScalarField& faceFlux = this->faceFlux_;

    // Check ENO stencil number storage
    if (stencilNo_.size() != sWeights.size())
    {
        selectStencil(vf);
    }

    forAll(sWeights, faceI)
    {
        bool success = false;

        scalar g = CDWeights[faceI]/(1.0 - CDWeights[faceI]);

        sWeights[faceI] = this->UENO2Weights(g, faceFlux[faceI], stencilNo_[faceI], success);
        if (success == false)
        {
            Pout<< "Invalid stencil index of internal face " << faceI << endl;
        }
    }

    surfaceScalarField::Boundary& bSWeights = sWeights.boundaryFieldRef();

    forAll(bSWeights, patchI)
    {
        fvsPatchField<scalar>& pSWeights = bSWeights[patchI];

        if (pSWeights.coupled())
        {
            const scalarField& pCDWeights = CDWeights.boundaryField()[patchI];

            const scalarField& pFaceFlux = faceFlux.boundaryField()[patchI];

            const labelField& pStencilNo = bStencilNo_[patchI];

            forAll(pSWeights, faceI)
            {
                bool success = false;

                scalar g = pCDWeights[faceI]/(1.0 - pCDWeights[faceI]);

                pSWeights[faceI] = this->UENO2Weights(g, pFaceFlux[faceI], pStencilNo[faceI], success);
                if (success == false)
                {
                    Pout<< "Invalid stencil index of patch " << pSWeights.patch().name() << " face " << faceI << endl;
                }
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
Foam::UENO2<Type>::correction
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
                "UENO2::correction(" + vf.name() + ')',
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

    const surfaceScalarField& faceFlux = this->faceFlux_;

    // Check ENO stencil number storage
    if (stencilNo_.size() != sfCorr.size())
    {
        selectStencil(vf);
    }

    const GeometricField<Type, fvsPatchField, surfaceMesh>& faceOwnDiff1 = this->faceOwnDiff1();
    const GeometricField<Type, fvsPatchField, surfaceMesh>& faceNeiDiff1 = this->faceNeiDiff1();

    forAll(sfCorr, faceI)
    {
        bool success = false;

        scalar g = CDWeights[faceI]/(1.0 - CDWeights[faceI]);

        sfCorr[faceI] =
            this->UENO2Corr
            (
                g,
                faceOwnDiff1[faceI],
                faceNeiDiff1[faceI],
                faceFlux[faceI],
                stencilNo_[faceI],
                success
            );
        if (success == false)
        {
            Pout<< "Invalid stencil index of internal face " << faceI << endl;
        }
    }

    typename GeometricField<Type, fvsPatchField, surfaceMesh>::
        Boundary& bSfCorr = sfCorr.boundaryFieldRef();

    forAll(bSfCorr, patchI)
    {
        fvsPatchField<Type>& pSfCorr = bSfCorr[patchI];

        if (pSfCorr.coupled())
        {
            const scalarField& pCDWeights = CDWeights.boundaryField()[patchI];

            const scalarField& pFaceFlux = faceFlux.boundaryField()[patchI];

            const labelField& pStencilNo = bStencilNo_[patchI];

            const Field<Type>& pFaceOwnDiff1 = faceOwnDiff1.boundaryField()[patchI];
            const Field<Type>& pFaceNeiDiff1 = faceNeiDiff1.boundaryField()[patchI];

            forAll(pSfCorr, faceI)
            {
                bool success = false;

                scalar g = pCDWeights[faceI]/(1.0 - pCDWeights[faceI]);

                pSfCorr[faceI] =
                    this->UENO2Corr
                    (
                        g,
                        pFaceOwnDiff1[faceI],
                        pFaceNeiDiff1[faceI],
                        pFaceFlux[faceI],
                        pStencilNo[faceI],
                        success
                    );
                if (success == false)
                {
                    Pout<< "Invalid stencil index of patch " << pSfCorr.patch().name() << " face " << faceI << endl;
                }
            }
        }
        else
        {
            forAll(pSfCorr, faceI)
            {
                pSfCorr[faceI] = Type(Zero);
            }
        }
    }

    return tsfCorr;
}


template<class Type>
Foam::scalar Foam::UENO2<Type>::UENO2Weights
(
    scalar g,
    scalar faceFlux,
    label stencilNo,
    bool& success
) const
{
    success = true;

    if (stencilNo == 0)
    {
        return  g/(1 + g);
    }
    else if (stencilNo == 1)
    {
        if (faceFlux >= 0)
        {
            scalar c10 = -g/(1 + g);
            scalar c11 = 1 - c10;

            return c10*this->cellL1Weights(g) + c11;
        }
        else
        {
            g = 1.0/g;
            scalar c10 = -g/(1 + g);
            scalar c11 = 1 - c10;

            return 1.0 - (c10*this->cellL1Weights(g) + c11);
        }
    }
    else
    {
        success = false;
        return 0.5;
    }

}


template<class Type>
Type Foam::UENO2<Type>::UENO2Corr
(
    scalar g,
    const Type& faceOwnDiff1,
    const Type& faceNeiDiff1,
    scalar faceFlux,
    label stencilNo,
    bool& success
) const
{
    success = true;

    if (stencilNo == 0)
    {
        return Type(Zero);
    }
    else if (stencilNo == 1)
    {
        if (faceFlux >= 0)
        {
            scalar c10 = -g/(1 + g);

            return c10*this->cellL1Corr(faceOwnDiff1, g);
        }
        else
        {
            g = 1.0/g;
            scalar c10 = -g/(1 + g);

            return c10*this->cellL1Corr(-faceNeiDiff1, g);
        }
    }
    else
    {
        success = false;
        return Type(Zero);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeFirstOrderDolincScheme(UENO2)
}

// ************************************************************************* //
