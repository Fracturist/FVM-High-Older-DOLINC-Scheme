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
    Foam::CSR4

Description
    4th-order central-stencil reconstruction scheme class

SourceFiles
    CSR4.C

Author
    Hao Guo

\*---------------------------------------------------------------------------*/

#ifndef CSR4_H
#define CSR4_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class CSR4 Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CSR4
:
    public surfaceInterpolationScheme<Type>
{
    // Private Data

        word gradSchemeName_;


    // Private Member Functions

        //- No copy assignment
        void operator=(const CSR4&) = delete;

        //- Calculate scheme weights using mesh grow rate
        inline scalar schemeWeights(scalar g) const
        {
            scalar c10
            (
                -pow5(g)
               /(1 + g + sqr(g))
               /(1 + g + sqr(g) + pow3(g))
            );
            scalar c11
            (
                sqr(g)
              *(2*pow3(g) + 2*sqr(g) + 2*g  + 1)
              /(1 + g + sqr(g))
              /(1 + g + sqr(g) + pow3(g))
            );
            scalar c13
            (
                -1
               /(1 + g + sqr(g))
               /(1 + g + sqr(g) + pow3(g))
            );

            return (1 - 1/g)*c10 + c11 + g*c13;
        }

        //- Calculate scheme correction part
        inline scalar schemeCorr
        (
            scalar g,
            const vector& gradO,
            const vector& gradN,
            const vector& delta
        ) const
        {
            scalar c10
            (
                -pow5(g)
               /(1 + g + sqr(g))
               /(1 + g + sqr(g) + pow3(g))
            );
            scalar c13
            (
                -1
               /(1 + g + sqr(g))
               /(1 + g + sqr(g) + pow3(g))
            );

            return (g*gradN*c13 - gradO*c10/g) & (2*delta);
        }


public:

    //- Runtime type information
    TypeName("CSR4");


    // Constructors

        //- Construct from mesh
        CSR4(const fvMesh& mesh)
        :
            surfaceInterpolationScheme<Type>(mesh),
            gradSchemeName_("grad")
        {}

        //- Construct from Istream
        CSR4(const fvMesh& mesh, Istream& schemeData)
        :
            surfaceInterpolationScheme<Type>(mesh),
            gradSchemeName_(schemeData)
        {}

        //- Construct from faceFlux and Istream
        CSR4
        (
            const fvMesh& mesh,
            const surfaceScalarField&,
            Istream& schemeData
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            gradSchemeName_(schemeData)
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
