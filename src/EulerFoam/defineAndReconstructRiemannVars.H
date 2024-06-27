surfaceScalarField rho_negL(interpolate(rho, negL));
surfaceScalarField rho_posR(interpolate(rho, posR));

surfaceVectorField rhoU_negL(interpolate(rhoU, negL, U.name()));
surfaceVectorField rhoU_posR(interpolate(rhoU, posR, U.name()));

surfaceScalarField p_negL(interpolate(p, negL, T.name()));
surfaceScalarField p_posR(interpolate(p, posR, T.name()));


surfaceVectorField U_negL("U_negL", rhoU_negL / rho_negL);
surfaceVectorField U_posR("U_posR", rhoU_posR / rho_posR);

surfaceScalarField phiv_negL("phiv_negL", U_negL & mesh.Sf());
phiv_negL.setOriented(false);
surfaceScalarField phiv_posR("phiv_posR", U_posR & mesh.Sf());
phiv_posR.setOriented(false);


volScalarField gamma(thermo.gamma());
surfaceScalarField gamma_negL(upwindProperty(gamma, negL));
surfaceScalarField gamma_posR(upwindProperty(gamma, posR));

surfaceScalarField rhoE_negL
(
    p_negL/(gamma_negL - 1.0)
  + 0.5*magSqr(rhoU_negL)/rho_negL
);
surfaceScalarField rhoE_posR
(
    p_posR/(gamma_posR - 1.0)
  + 0.5*magSqr(rhoU_posR)/rho_posR
);

volScalarField a("a", sqrt(gamma*p/rho));
surfaceScalarField aSf_negL
(
    "aSf_negL",
    sqrt(gamma_negL*p_negL/rho_negL)*mesh.magSf()
);
surfaceScalarField aSf_posR
(
    "aSf_posR",
    sqrt(gamma_posR*p_posR/rho_posR)*mesh.magSf()
);
