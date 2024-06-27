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
    Foam::USR3

Description
    One cell upwind biased 3rd-order fixed stencil reconstruction scheme class

SourceFiles
    USR3.C

Author
    Hao Guo

\*---------------------------------------------------------------------------*/

#ifndef USR3_H
#define USR3_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class USR3 Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class USR3
:
    public surfaceInterpolationScheme<Type>
{
    // Private Data

        //- Reference to the surface flux used to choose upwind direction
        const surfaceScalarField& faceFlux_;

        //- Gauss gradient scheme name
        word gradSchemeName_;


    // Private Member Functions

        //- No copy assignment
        void operator=(const USR3&) = delete;

        //- Calculate upwind stencil scheme weights using mesh grow rate
        inline scalar schemeWeights(scalar g) const
        {
            scalar c10
            (
                -pow3(g)
               /(1 + g)
               /(1 + g + sqr(g))
            );
            scalar c11
            (
                g
               *(1 + 2*g + 2*sqr(g))
               /(1 + g)
               /(1 + g + sqr(g))
            );

            return (1 - 1/g)*c10 + c11;
        }

        //- Calculate upwind stencil scheme correction part
        inline scalar schemeCorr
        (
            scalar g,
            const vector& gradP,
            const vector& delta
        ) const
        {
            scalar c10
            (
                -pow3(g)
               /(1 + g)
               /(1 + g + sqr(g))
            );

            return (-gradP*c10/g) & (2*delta);
        }


public:

    //- Runtime type information
    TypeName("USR3");


    // Constructors

        //- Construct from mesh and faceFlux
        USR3
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_(faceFlux),
            gradSchemeName_("grad")
        {}

        //- Construct from mesh and Istream
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        USR3(const fvMesh& mesh, Istream& is)
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_(mesh.lookupObject<surfaceScalarField>(word(is))),
            gradSchemeName_(is)
        {}

        //- Construct from faceFlux and Istream
        USR3
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_(faceFlux),
            gradSchemeName_(is)
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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
