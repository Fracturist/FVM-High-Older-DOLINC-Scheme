#include "upwind.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Interpolate field vf according to direction dir
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const surfaceScalarField& dir,
    const word& reconFieldName = word::null
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsf
    (
        fvc::interpolate
        (
            vf,
            dir,
            "reconstruct("
          + (reconFieldName != word::null ? reconFieldName : vf.name())
          + ')'
        )
    );

    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf.ref();

    sf.rename(vf.name() + '_' + dir.name());

    return tsf;
}

//- Interpolate thermophysical property field according to direction dir
//  using upwind scheme (cell value)
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> upwindProperty
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const surfaceScalarField& dir
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsf
    (
        upwind<Type>(dir.mesh(), dir).interpolate(vf)
    );

    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf.ref();

    sf.rename(vf.name() + '_' + dir.name());

    return tsf;
}

}
