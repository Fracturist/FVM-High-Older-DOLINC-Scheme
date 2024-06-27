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
    Foam::UENO3

Group
    grpFvSurfaceInterpolationSchemes

Description
    3rd order ENO scheme

SourceFiles
    UENO3.C

Author
    Hao Guo

\*---------------------------------------------------------------------------*/

#ifndef UENO3_H
#define UENO3_H

#include "secondOrderDolincScheme.H"
#include "volFields.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class UENO3 Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class UENO3
:
    public secondOrderDolincScheme<Type>
{
    // Private Data

        //- Limit factor used for stencil smoothness evaluation
        scalar b_;

        //- Index of stencil selected by ENO method
        mutable labelField stencilNo_;
        mutable List<labelField> bStencilNo_;


    // Private Member Functions

        //- No copy assignment
        void operator=(const UENO3&) = delete;

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

        //- Calculate smoothness of S_3( I_{i-2}, I_{i-1}, I_{i} )
        inline scalar smoothness3L2
        (
            scalar g,
            const Type& vfP,
            const Type& vfD,
            const Type& diffP1,
            const Type& diffP2,
            const Type& diffD1
        ) const
        {
            return
            mag
            (
                sqr(g)
               *(
                    vfP
                  - (1.0 + g)*this->vfL1(vfP, vfD, diffP1, g)
                  + g*this->vfL2(vfP, vfD, diffP1, diffP2, diffD1, g)
                )
            );
        }

        //- Calculate smoothness of S_3( I_{i-1}, I_{i}, I_{i+1} )
        inline scalar smoothness3L1
        (
            scalar g,
            const Type& vfP,
            const Type& vfD,
            const Type& diffP1,
            const Type& diffD1
        ) const
        {
            return
            mag
            (
                vfD
              - (1.0 + g)*vfP
              + g*this->vfL1(vfP, vfD, diffP1, g)
            );
        }

        //- Calculate smoothness of S_3( I_{i}, I_{i+1}, I_{i+2} )
        inline scalar smoothness3L0
        (
            scalar g,
            const Type& vfP,
            const Type& vfD,
            const Type& diffD1
        ) const
        {
            return
            mag
            (
                (
                    this->vfR2(vfP, vfD, diffD1, g)
                  - (1.0 + g)*vfD
                  + g*vfP
                )
               /sqr(g)
            );
        }


public:

    //- Runtime type information
    TypeName("UENO3");


    // Constructors

        //- Construct from Istream
        UENO3(const fvMesh& mesh, Istream& is)
        :
            secondOrderDolincScheme<Type>(mesh, is),
            b_(readScalar(is))
        {}

        //- Construct from faceFlux and Istream
        UENO3
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            secondOrderDolincScheme<Type>(mesh, faceFlux, is),
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
        scalar UENO3Weights
        (
            scalar g,
            scalar faceFlux,
            label stencilNo,
            bool& success
        ) const;

        //- Return scheme correction of upwind ENO2
        Type UENO3Corr
        (
            scalar g,
            const Type& faceOwnDiff1,
            const Type& faceOwnDiff2,
            const Type& faceNeiDiff1,
            const Type& faceNeiDiff2,
            scalar faceFlux,
            label stencilNo,
            bool& success
        ) const;

        //- Return G(n)
        inline scalar G
        (
            scalar g,
            label n
        ) const
        {
            scalar res = 1.0;

            while(n > 1)
            {
                res = res*g + 1.0;
                n -= 1;
            }

            return res;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
