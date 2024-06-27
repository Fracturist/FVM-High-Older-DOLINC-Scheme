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
    Foam::secondOrderDolincScheme

Group
    grpFvLimitedSurfaceInterpolationSchemes

Description
    Abstract base class for limited surface interpolation schemes.

SourceFiles
    secondOrderDolincScheme.C

Author
    Hao Guo

\*---------------------------------------------------------------------------*/

#ifndef secondOrderDolincScheme_H
#define secondOrderDolincScheme_H

#include "surfaceInterpolationScheme.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class secondOrderDolincScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class secondOrderDolincScheme
:
    public surfaceInterpolationScheme<Type>
{
    // Private Member Functions

        //- No copy construct
        secondOrderDolincScheme
        (
            const secondOrderDolincScheme&
        ) = delete;

        //- No copy assignment
        void operator=(const secondOrderDolincScheme&) = delete;


protected:

    // Protected data

        const surfaceScalarField& faceFlux_;

        //- Gauss gradient scheme name
        word gradSchemeName_;

        //- Inner product of owner or neighbour cell gradient and face 2 delta
        mutable autoPtr<GeometricField<Type, fvsPatchField, surfaceMesh>> faceOwnDiff1_;
        mutable autoPtr<GeometricField<Type, fvsPatchField, surfaceMesh>> faceNeiDiff1_;
        mutable autoPtr<GeometricField<Type, fvsPatchField, surfaceMesh>> faceOwnDiff2_;
        mutable autoPtr<GeometricField<Type, fvsPatchField, surfaceMesh>> faceNeiDiff2_;


    // Protected Member Functions

        //- Calculate and store differentiation of several orders
        void calcDiff
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;


public:

    //- Runtime type information
    TypeName("secondOrderDolincScheme");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            secondOrderDolincScheme,
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
            secondOrderDolincScheme,
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
        secondOrderDolincScheme
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
        secondOrderDolincScheme
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
        static tmp<secondOrderDolincScheme<Type>> New
        (
            const fvMesh& mesh,
            Istream& schemeData
        );

        //- Return new tmp interpolation scheme
        static tmp<secondOrderDolincScheme<Type>> New
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        );


    //- Destructor
    virtual ~secondOrderDolincScheme();


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
             || !faceOwnDiff2_.valid()
             || !faceNeiDiff1_.valid()
             || !faceNeiDiff2_.valid()
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

        //- Return owner side 2nd order difference field
        inline const GeometricField<Type, fvsPatchField, surfaceMesh>&
        faceOwnDiff2() const
        {
            return faceOwnDiff2_();

        }

        //- Return neighbour side 1st order difference field
        inline const GeometricField<Type, fvsPatchField, surfaceMesh>&
        faceNeiDiff1() const
        {
            return faceNeiDiff1_();
        }

        //- Return neighbour side 2nd order difference field
        inline const GeometricField<Type, fvsPatchField, surfaceMesh>&
        faceNeiDiff2() const
        {
            return faceNeiDiff2_();
        }


        // Dolinc accesses

            //- Return vf value of cell i-2
            inline Type vfL2
            (
                const Type& vfP,
                const Type& vfD,
                const Type& diffP1,
                const Type& diffP2,
                const Type& diffD1,
                scalar g = 1.0
            ) const
            {
                return
                (
                    ( 1.0 - (1.0 - 1.0/g)/g )*vfP
                  + (1.0 - 1.0/g)/g*vfD
                  - (1.0 - 1.0/sqr(g))/g*diffP1
                  + (diffP2 - diffD1)/pow3(g)
                );
            }

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

            //- Return vf value of cell i+3
            inline Type vfR3
            (
                const Type& vfP,
                const Type& vfD,
                const Type& diffP1,
                const Type& diffD1,
                const Type& diffD2,
                scalar g = 1.0
            ) const
            {
                return
                (
                    (1.0 - g)*g*vfP
                  + ( 1.0 - (1.0 - g)*g )*vfD
                  + (1.0 - sqr(g))*g*diffD1
                  + pow3(g)*(diffD2 + diffP1)
                );
            }

            //- Return weights of cell i-2 DOLINC
            inline scalar cellL2Weights(scalar g = 1.0) const
            {
                return 1.0 - (1.0 - 1.0/g)/g;
            }

            //- Return correction of cell i-2 DOLINC
            inline Type cellL2Corr
            (
                const Type& diffP1,
                const Type& diffP2,
                const Type& diffD1,
                scalar g = 1.0
            ) const
            {
                return -(1.0 - 1.0/sqr(g))/g*diffP1 + (diffP2 - diffD1)/pow3(g);
            }

            //- Return correction of cell i-2 DOLINC using cmpt
            inline scalar cellL2CorrCmpt
            (
                scalar diffP1Cmpt,
                scalar diffP2Cmpt,
                scalar diffD1Cmpt,
                scalar g = 1.0
            ) const
            {
                return -(1.0 - 1.0/sqr(g))/g*diffP1Cmpt + (diffP2Cmpt - diffD1Cmpt)/pow3(g);
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

            //- Return weights of cell i+3 DOLINC
            inline scalar cellR3Weights(scalar g = 1.0) const
            {
                return g*(1.0 - g);
            }

            //- Return correction of cell i+3 DOLINC
            inline Type cellR3Corr
            (
                const Type& diffP1,
                const Type& diffD1,
                const Type& diffD2,
                scalar g = 1.0
            ) const
            {
                return (1.0 - sqr(g))*g*diffD1 + pow3(g)*(diffD2 + diffP1);
            }

            //- Return correction of cell i+3 DOLINC using cmpt
            inline scalar cellR3CorrCmpt
            (
                scalar diffP1Cmpt,
                scalar diffD1Cmpt,
                scalar diffD2Cmpt,
                scalar g = 1.0
            ) const
            {
                return (1.0 - sqr(g))*g*diffD1Cmpt + pow3(g)*(diffD2Cmpt + diffP1Cmpt);

            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeSecondOrderDolincTypeScheme(SS, Type)                    \
                                                                               \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                              \
                                                                               \
surfaceInterpolationScheme<Type>::addMeshConstructorToTable<SS<Type>>          \
    add##SS##Type##MeshConstructorToTable_;                                    \
                                                                               \
surfaceInterpolationScheme<Type>::addMeshFluxConstructorToTable<SS<Type>>      \
    add##SS##Type##MeshFluxConstructorToTable_;                                \
                                                                               \
secondOrderDolincScheme<Type>::addMeshConstructorToTable<SS<Type>>   \
    add##SS##Type##MeshConstructorToLimitedTable_;                             \
                                                                               \
secondOrderDolincScheme<Type>::                                      \
    addMeshFluxConstructorToTable<SS<Type>>                                    \
    add##SS##Type##MeshFluxConstructorToLimitedTable_;

#define makeSecondOrderDolincScheme(SS)                              \
                                                                               \
makeSecondOrderDolincTypeScheme(SS, scalar)                          \
makeSecondOrderDolincTypeScheme(SS, vector)                          \
makeSecondOrderDolincTypeScheme(SS, sphericalTensor)                 \
makeSecondOrderDolincTypeScheme(SS, symmTensor)                      \
makeSecondOrderDolincTypeScheme(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "secondOrderDolincScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
