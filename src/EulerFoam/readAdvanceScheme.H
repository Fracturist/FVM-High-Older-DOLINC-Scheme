autoPtr<volScalarField> aRhoc;
autoPtr<volVectorField> aRhoUc;
autoPtr<volScalarField> aRhoEc;

word advanceScheme("ForwardEuler");

if (mesh.schemesDict().readIfPresent("advanceScheme", advanceScheme))
{
    if
    (
        (advanceScheme == "ForwardEuler")
    )
    {
        Info<< "\nadvanceScheme: " << advanceScheme << endl;
    }
    else if
    (
        (advanceScheme == "TVDRungeKutta3")
    )
    {
        aRhoc = autoPtr<volScalarField>::New(rho);
        aRhoUc = autoPtr<volVectorField>::New(rhoU);
        aRhoEc = autoPtr<volScalarField>::New(rhoE);

        Info<< "\nadvanceScheme: " << advanceScheme << endl;
    }
    else
    {
        FatalErrorInFunction
            << "advanceScheme: " << advanceScheme
            << " is not a valid choice. "
            << "Options are: ForwardEuler, TVDRungeKutta3 and RungeKutta4"
            << abort(FatalError);
    }
}
