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
    Foam::firstOrderDolincScheme

Group
    grpFvLimitedSurfaceInterpolationSchemes

Description
    Abstract base class for limited surface interpolation schemes.

SourceFiles
    firstOrderDolincScheme.C

Author
    Hao Guo

\*---------------------------------------------------------------------------*/

#ifndef firstOrderDolincScheme_H
#define firstOrderDolincScheme_H

#include "surfaceInterpolationScheme.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class firstOrderDolincScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class firstOrderDolincScheme
:
    public surfaceInterpolationScheme<Type>
{
    // Private Member Functions

        //- No copy construct
        firstOrderDolincScheme
        (
            const firstOrderDolincScheme&
        ) = delete;

        //- No copy assignment
        void operator=(const firstOrderDolincScheme&) = delete;


protected:

    // Protected data

        const surfaceScalarField& faceFlux_;

        //- Gauss gradient scheme name
        word gradSchemeName_;

        //- Inner product of owner or neighbour cell gradient and face 2 delta
        mutable autoPtr<GeometricField<Type, fvsPatchField, surfaceMesh>> faceOwnDiff1_;
        mutable autoPtr<GeometricField<Type, fvsPatchField, surfaceMesh>> faceNeiDiff1_;


    // Protected Member Functions

        //- Calculate and store differentiation of several orders
        void calcDiff
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;


public:

    //- Runtime type information
    TypeName("firstOrderDolincScheme");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            firstOrderDolincScheme,
            Mesh,
            (
                const fvMesh& mesh,
                Istream& schemeData
            ),
            (mesh, schemeData)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            firstOrderDolincScheme,
            MeshFlux,
            (
                const fvMesh& mesh,
                const surfaceScalarField& faceFlux,
                Istream& schemeData
            ),
            (mesh, faceFlux, schemeData)
        );


    // Constructors

        //- Construct from mesh and Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        firstOrderDolincScheme
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_
            (
                mesh.lookupObject<surfaceScalarField>
                (
                    word(is)
                )
            ),
            gradSchemeName_(is)
        {}


        //- Construct from mesh, faceFlux and Istream
        firstOrderDolincScheme
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


    // Selectors

        //- Return new tmp interpolation scheme
        static tmp<firstOrderDolincScheme<Type>> New
        (
            const fvMesh& mesh,
            Istream& schemeData
        );

        //- Return new tmp interpolation scheme
        static tmp<firstOrderDolincScheme<Type>> New
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        );


    //- Destructor
    virtual ~firstOrderDolincScheme();


    // Member Functions

        //- Check differentiation storage
        inline void checkDiff
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            if
            (
                !faceOwnDiff1_.valid()
             || !faceNeiDiff1_.valid()
            )
            {
                calcDiff(vf);
            }
        }


        //- Return owner side 1st order difference field
        inline const GeometricField<Type, fvsPatchField, surfaceMesh>&
        faceOwnDiff1() const
        {
            return faceOwnDiff1_();

        }

        //- Return neighbour side 1st order difference field
        inline const GeometricField<Type, fvsPatchField, surfaceMesh>&
        faceNeiDiff1() const
        {
            return faceNeiDiff1_();
        }


        // Dolinc accesses

            //- Return vf value of cell i-1
            inline Type vfL1
            (
                const Type& vfP,
                const Type& vfD,
                const Type& diffP1,
                scalar g = 1.0
            ) const
            {
                return 1.0/g*vfD + (g-1.0)/g*vfP - 1.0/g*diffP1;
            }

            //- Return vf value of cell i+2
            inline Type vfR2
            (
                const Type& vfP,
                const Type& vfD,
                const Type& diffD1,
                scalar g = 1.0
            ) const
            {
                return g*vfP + (1.0-g)*vfD + g*diffD1;
            }

            //- Return weights of cell i-1 DOLINC
            inline scalar cellL1Weights(scalar g = 1.0) const
            {
                return 1.0 - 1.0/g;
            }

            //- Return correction of cell i-1 DOLINC
            inline Type cellL1Corr
            (
                const Type& diffP1,
                scalar g = 1.0
            ) const
            {
                return -1.0/g*diffP1;
            }

            //- Return correction of cell i-1 DOLINC using cmpt
            inline scalar cellL1CorrCmpt
            (
                scalar diffP1Cmpt,
                scalar g = 1.0
            ) const
            {
                return -1.0/g*diffP1Cmpt;
            }

            //- Return weights of cell i+2 DOLINC
            inline scalar cellR2Weights(scalar g = 1.0) const
            {
                return g;
            }

            //- Return correction of cell i+2 DOLINC
            inline Type cellR2Corr
            (
                const Type& diffD1,
                scalar g = 1.0
            ) const
            {
                return g*diffD1;
            }

            //- Return correction of cell i+2 DOLINC using cmpt
            inline scalar cellR2CorrCmpt
            (
                scalar diffD1Cmpt,
                scalar g = 1.0
            ) const
            {
                return g*diffD1Cmpt;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeFirstOrderDolincTypeScheme(SS, Type)                    \
                                                                               \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                              \
                                                                               \
surfaceInterpolationScheme<Type>::addMeshConstructorToTable<SS<Type>>          \
    add##SS##Type##MeshConstructorToTable_;                                    \
                                                                               \
surfaceInterpolationScheme<Type>::addMeshFluxConstructorToTable<SS<Type>>      \
    add##SS##Type##MeshFluxConstructorToTable_;                                \
                                                                               \
firstOrderDolincScheme<Type>::addMeshConstructorToTable<SS<Type>>   \
    add##SS##Type##MeshConstructorToLimitedTable_;                             \
                                                                               \
firstOrderDolincScheme<Type>::                                      \
    addMeshFluxConstructorToTable<SS<Type>>                                    \
    add##SS##Type##MeshFluxConstructorToLimitedTable_;

#define makeFirstOrderDolincScheme(SS)                              \
                                                                               \
makeFirstOrderDolincTypeScheme(SS, scalar)                          \
makeFirstOrderDolincTypeScheme(SS, vector)                          \
makeFirstOrderDolincTypeScheme(SS, sphericalTensor)                 \
makeFirstOrderDolincTypeScheme(SS, symmTensor)                      \
makeFirstOrderDolincTypeScheme(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "firstOrderDolincScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
