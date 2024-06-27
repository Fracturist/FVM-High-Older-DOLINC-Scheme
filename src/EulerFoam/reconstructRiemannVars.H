rho_negL = interpolate(rho, negL);
rho_posR = interpolate(rho, posR);

rhoU_negL = interpolate(rhoU, negL, U.name());
rhoU_posR = interpolate(rhoU, posR, U.name());

p_negL = interpolate(p, negL, T.name());
p_posR = interpolate(p, posR, T.name());


U_negL = rhoU_negL / rho_negL;
U_posR = rhoU_posR / rho_posR;

phiv_negL = U_negL & mesh.Sf();
phiv_negL.setOriented(false);
phiv_posR = U_posR & mesh.Sf();
phiv_posR.setOriented(false);


gamma = thermo.gamma();
gamma_negL = upwindProperty(gamma, negL);
gamma_posR = upwindProperty(gamma, posR);

rhoE_negL =
    p_negL/(gamma_negL - 1.0)
  + 0.5*magSqr(rhoU_negL)/rho_negL;
rhoE_posR =
    p_posR/(gamma_posR - 1.0)
  + 0.5*magSqr(rhoU_posR)/rho_posR;

a = sqrt(gamma*p/rho);
aSf_negL =
    sqrt(gamma_negL*p_negL/rho_negL)*mesh.magSf();
aSf_posR =
    sqrt(gamma_posR*p_posR/rho_posR)*mesh.magSf();
