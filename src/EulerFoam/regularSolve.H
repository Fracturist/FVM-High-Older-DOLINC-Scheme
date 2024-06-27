// --- Solve density
solve(fvm::ddt(rho) + fvc::div(fluxRho));

// --- Solve momentum
solve(fvm::ddt(rhoU) + fvc::div(fluxRhoU));

// --- Solve energy
solve(fvm::ddt(rhoE) + fvc::div(fluxRhoE));


U.ref() = rhoU()/rho();
U.correctBoundaryConditions();

e.ref() = rhoE()/rho() - 0.5*magSqr(U());
e.correctBoundaryConditions();
thermo.correct();

p.ref() = rho()/psi();
p.correctBoundaryConditions();

rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();
rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();
rhoE.boundaryFieldRef() ==
    rho.boundaryField()*
    (
        e.boundaryField() + 0.5*magSqr(U.boundaryField())
    );
