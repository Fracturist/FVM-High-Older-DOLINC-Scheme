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

Class
    Foam::UENO2

Group
    grpFvSurfaceInterpolationSchemes

Description
    2nd order ENO scheme

SourceFiles
    UENO2.C

Author
    Hao Guo

\*---------------------------------------------------------------------------*/

#ifndef UENO2_H
#define UENO2_H

#include "firstOrderDolincScheme.H"
#include "volFields.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class UENO2 Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class UENO2
:
    public firstOrderDolincScheme<Type>
{
    // Private Data

        //- Limit factor used for stencil smoothness evaluation
        scalar b_;

        //- Index of stencil selected by ENO method
        mutable labelField stencilNo_;
        mutable List<labelField> bStencilNo_;


    // Private Member Functions

        //- No copy assignment
        void operator=(const UENO2&) = delete;

        //- Calculate stencil smoothness and set the selected index
        void selectStencil
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;

        //- Calculate smoothness of S_2( I_{i-1}, I_{i} )
        inline scalar smoothness2L1
        (
            scalar g,
            const Type& vfP,
            const Type& vfD,
            const Type& diffP1
        ) const
        {
            return mag( g*(vfP - this->vfL1(vfP, vfD, diffP1, g)) );
        }

        //- Calculate smoothness of S_2( I_{i}, I_{i+1} )
        inline scalar smoothness2L0
        (
            const Type& vfP,
            const Type& vfD
        ) const
        {
            return mag( vfD - vfP );
        }


public:

    //- Runtime type information
    TypeName("UENO2");


    // Constructors

        //- Construct from Istream
        UENO2(const fvMesh& mesh, Istream& is)
        :
            firstOrderDolincScheme<Type>(mesh, is),
            b_(readScalar(is))
        {}

        //- Construct from faceFlux and Istream
        UENO2
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            firstOrderDolincScheme<Type>(mesh, faceFlux, is),
            b_(readScalar(is))
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Return scheme weights of upwind ENO2
        scalar UENO2Weights
        (
            scalar g,
            scalar faceFlux,
            label stencilNo,
            bool& success
        ) const;

        //- Return scheme correction of upwind ENO2
        Type UENO2Corr
        (
            scalar g,
            const Type& faceOwnDiff1,
            const Type& faceNeiDiff1,
            scalar faceFlux,
            label stencilNo,
            bool& success
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
