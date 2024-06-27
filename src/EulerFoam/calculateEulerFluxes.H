if (fluxScheme == "Tadmor")
{
    const scalar w_ab = 0.5;

    surfaceScalarField lambda_ab
    (
        "lambda_ab",
        max
        (
            max
            (
                mag(phiv_negL + aSf_negL),
                mag(phiv_posR + aSf_posR)
            ),
            max
            (
                mag(phiv_negL - aSf_negL),
                mag(phiv_posR - aSf_posR)
            )
        )
    );

    fluxRho =
        w_ab*(phiv_negL + lambda_ab)*rho_negL
      + w_ab*(phiv_posR - lambda_ab)*rho_posR;

    fluxRhoU =
        w_ab*(phiv_negL + lambda_ab)*rhoU_negL
      + w_ab*(phiv_posR - lambda_ab)*rhoU_posR;
    fluxRhoU.setOriented(true);
    fluxRhoU += w_ab*(p_negL + p_posR)*mesh.Sf();

    fluxRhoE =
        w_ab*(phiv_negL*(rhoE_negL + p_negL) + lambda_ab*rhoE_negL)
      + w_ab*(phiv_posR*(rhoE_posR + p_posR) - lambda_ab*rhoE_posR);
}
else if (fluxScheme == "Kurganov")
{
    const dimensionedScalar v_zero(dimVolume/dimTime, Zero);

    surfaceScalarField aMin
    (
        "aMin",
        min(min(phiv_negL - aSf_negL, phiv_posR - aSf_posR), -v_zero)
    );

    surfaceScalarField aMax
    (
        "aMax",
        max(max(phiv_negL + aSf_negL, phiv_posR + aSf_posR), v_zero)
    );

    surfaceScalarField w_a
    (
        "w_a",
        aMax/(aMax - aMin)
    );

    surfaceScalarField w_b
    (
        "w_b",
        1.0 - w_a
    );

    fluxRho =
        w_a*(phiv_negL - aMin)*rho_negL
      + w_b*(phiv_posR - aMax)*rho_posR;

    fluxRhoU =
        w_a*(phiv_negL - aMin)*rhoU_negL
      + w_b*(phiv_posR - aMax)*rhoU_posR;
    fluxRhoU.setOriented(true);
    fluxRhoU += (w_a*p_negL + w_b*p_posR)*mesh.Sf();

    fluxRhoE =
        w_a*(phiv_negL*(rhoE_negL + p_negL) - aMin*rhoE_negL)
      + w_b*(phiv_posR*(rhoE_posR + p_posR) - aMax*rhoE_posR);
}
