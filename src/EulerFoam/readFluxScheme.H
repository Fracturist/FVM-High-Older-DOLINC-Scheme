word fluxScheme("Kurganov");
if (mesh.schemesDict().readIfPresent("fluxScheme", fluxScheme))
{
    if
    (
        (fluxScheme == "Tadmor")
     || (fluxScheme == "Kurganov")
    )
    {
        Info<< "\nfluxScheme: " << fluxScheme << endl;
    }
    else
    {
        FatalErrorInFunction
            << "fluxScheme: " << fluxScheme
            << " is not a valid choice. "
            << "Options are: Lax, Friedrichs, Tadmor, Kurganov and Roe"
            << abort(FatalError);
    }
}
