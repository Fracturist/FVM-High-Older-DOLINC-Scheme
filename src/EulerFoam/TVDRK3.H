{
    volScalarField& rhoc = aRhoc();
    volVectorField& rhoUc = aRhoUc();
    volScalarField& rhoEc = aRhoEc();
    rhoc = rho.oldTime();
    rhoUc = rhoU.oldTime();
    rhoEc = rhoE.oldTime();


    // * * * * * * * * * * * * * * U(1)  * * * * * * * * * * * * * * //

    Info<< ">>> 1st solution for U(1) >>>" << endl;

    // Un + deltaT*L(Un)

    #include "calculateEulerFluxes.H"
    #include "regularSolve.H"


    // * * * * * * * * * * * * * * U(2)  * * * * * * * * * * * * * * //

    Info<< ">>> 2nd solution for U(2) >>>" << endl;

    // 3/4*Un + 1/4*(U(1) + deltaT*L(U(1)))

    #include "reconstructRiemannVars.H"
    #include "calculateEulerFluxes.H"
    #include "fvmDdtSolve.H"

    rho = 3.0/4.0*rhoc + 1.0/4.0*rho;
    rhoU = 3.0/4.0*rhoUc + 1.0/4.0*rhoU;
    rhoE = 3.0/4.0*rhoEc + 1.0/4.0*rhoE;


    // * * * * * * * * * * * * * * Un+1  * * * * * * * * * * * * * * //

    Info<< ">>> Final solution for Un+1 >>>" << endl;

    // 1/3*Un + 2/3*(U(2) + deltaT*L(U(2)))

    #include "reconstructRiemannVars.H"
    #include "calculateEulerFluxes.H"
    #include "fvmDdtSolve.H"

    rho = 1.0/3.0*rhoc + 2.0/3.0*rho;
    rhoU = 1.0/3.0*rhoUc + 2.0/3.0*rhoU;
    rhoE = 1.0/3.0*rhoEc + 2.0/3.0*rhoE;
}
