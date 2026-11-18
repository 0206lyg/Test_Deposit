// Reconstruct one Cartesian normal component from the staggered fluxes of a
// rectangular cell. On this mesh each wall sample is at the midpoint of its
// tangential face bin. Integration over the bin recovers the interpolated flux.
#ifndef planarNormalVelocity_H
#define planarNormalVelocity_H

inline Foam::scalar planarNormalVelocity
(
    const Foam::fvMesh& mesh,
    const Foam::surfaceScalarField& phi,
    const Foam::label celli,
    const Foam::point& sample,
    const Foam::direction axis
)
{
    using namespace Foam;
    scalar coordinate[2] = {0, 0}, velocity[2] = {0, 0};
    bool found[2] = {false, false};
    const cell& faces = mesh.cells()[celli];
    forAll(faces, i)
    {
        const label face = faces[i];
        const vector& sf = mesh.faceAreas()[face];
        if (mag(sf[axis]) < 0.99*mag(sf)) continue;
        const scalar facePosition = mesh.faceCentres()[face][axis];
        const label side = facePosition > mesh.C()[celli][axis] ? 1 : 0;
        if (found[side])
            FatalErrorInFunction << "Wall flux interpolation needs a rectangular cell: "
                                 << celli << exit(FatalError);
        found[side] = true;
        coordinate[side] = facePosition;
        scalar flux;
        if (mesh.isInternalFace(face))
            flux = phi[face];
        else
        {
            const label patch = mesh.boundaryMesh().whichPatch(face);
            flux = phi.boundaryField()[patch][face - mesh.boundaryMesh()[patch].start()];
        }
        // phi and Sf share the same mesh orientation on internal, physical
        // and processor faces. Their ratio is signed along the global axis.
        velocity[side] = flux/sf[axis];
    }
    if (!found[0] || !found[1] || coordinate[1] <= coordinate[0])
        FatalErrorInFunction << "Missing bracketing faces for wall flux interpolation in cell "
                             << celli << exit(FatalError);
    const scalar fraction = (sample[axis]-coordinate[0])/(coordinate[1]-coordinate[0]);
    if (fraction < -1e-8 || fraction > 1+1e-8)
        FatalErrorInFunction << "Offset sample " << sample << " outside flux cell " << celli
                             << exit(FatalError);
    return velocity[0] + fraction*(velocity[1] - velocity[0]);
}
#endif
