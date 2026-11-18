// Cartesian outlet: traction_t = eta*(d_n U_t + grad_t U_n) + tau_nt.
// Incompressibility gives d_n U_n = -div_t U_t; normal traction fixes p.
{
    const fvPatch& patch = mesh.boundary()[outletID];
    const vectorField normals(patch.nf());
    const tmp<volTensorField> grad = fvc::grad(U);
    const tensorField& g = grad().boundaryField()[outletID];
    const tmp<volSymmTensorField> polymer = constEq.tau();
    const symmTensorField& tau = polymer().boundaryField()[outletID];
    fixedGradientFvPatchVectorField& outU =
        refCast<fixedGradientFvPatchVectorField>(U.boundaryFieldRef()[outletID]);
    fvPatchScalarField& outP = p.boundaryFieldRef()[outletID];
    // The constitutive model stores the updated apparent viscosity for Cross.
    const volScalarField* crossEta = depModel == "cross"
        ? &mesh.lookupObject<volScalarField>("eta") : nullptr;
    scalarField pressure(outP);
    forAll(normals, i)
    {
        const vector& n = normals[i];
        const scalar eta = crossEta ? crossEta->boundaryField()[outletID][i] : depEta;
        const scalar normalGrad = (n & (g[i] & n)) - tr(g[i]);
        const vector a = (g[i] & n) + (tau[i] & n)/eta;
        const vector gradient = normalGrad*n - (a - n*(n & a));
        outU.gradient()[i] += outletRelax*(gradient - outU.gradient()[i]);
        const scalar pn = (2*eta*normalGrad + (n & (tau[i] & n)))/constEq.rho().value();
        pressure[i] += outletRelax*(pn-pressure[i]);
    }
    outP == pressure;
    outU.evaluate();
}
