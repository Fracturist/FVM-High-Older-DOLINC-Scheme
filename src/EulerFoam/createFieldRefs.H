volScalarField& p = thermo.p();
const volScalarField& T = thermo.T();
const volScalarField& psi = thermo.psi();


bool debugF = args.found("debugFields");

if (debugF)
{
    e.write();
    psi.write();

    rhoU.write();
    rhoE.write();
}
