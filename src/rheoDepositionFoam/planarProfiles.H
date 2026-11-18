// Project-local diagnostics on orthogonal rectangular, one-cell-deep meshes.
namespace
{
using namespace Foam;
typedef FixedList<scalar, 5> WallRow; // wall ID, s, x, y, ds
const word wallPatches[] = {"wallsUpstream", "wallsContraction", "wallsDownstream"};
const char* wallNames[] = {"upstream", "step", "downstream"};

template<class Type>
scalar equationResidual(const fvMesh& mesh, const word& field)
{
    return cmptMax(Residuals<Type>::field(mesh, field)[0].initialResidual());
}

template<class Row>
List<List<Row>> collectRows(const UList<Row>& rows)
{
    List<List<Row>> blocks(Pstream::nProcs());
    blocks[Pstream::myProcNo()] = rows;
    Pstream::gatherList(blocks);
    return blocks;
}

void writeCells
(
    const fvMesh& mesh, const volVectorField& U, const surfaceScalarField& phi,
    const scalar depth, const label iteration, const fileName& directory
)
{
    typedef FixedList<scalar, 7> Row;
    List<Row> local(mesh.nCells());
    forAll(local, i)
    {
        const point& centre = mesh.C()[i];
        const cell& faces = mesh.cells()[i];
        scalar dx = 0, dy = 0, leftFlux = 0;
        forAll(faces, j)
        {
            const label f = faces[j];
            const vector& sf = mesh.faceAreas()[f];
            const point& cf = mesh.faceCentres()[f];
            dx = max(dx, 2*mag(cf.x()-centre.x()));
            dy = max(dy, 2*mag(cf.y()-centre.y()));
            if (mag(sf.x()) > 0.9*mag(sf) && cf.x() < centre.x())
            {
                if (mesh.isInternalFace(f))
                    leftFlux = phi[f];
                else
                {
                    const label patch = mesh.boundaryMesh().whichPatch(f);
                    leftFlux = phi.boundaryField()[patch][f-mesh.boundaryMesh()[patch].start()];
                }
                leftFlux *= sign(sf.x())/depth;
            }
        }
        Row& r = local[i];
        r[0]=centre.x(); r[1]=centre.y(); r[2]=dx; r[3]=dy;
        r[4]=U[i].x(); r[5]=U[i].y(); r[6]=leftFlux;
    }
    const List<List<Row>> cells = collectRows<Row>(local);
    if (Pstream::master())
    {
        std::ofstream out((directory/"flow_cells.csv").c_str());
        out.precision(16);
        out << "iteration,x_m,y_m,dx_m,dy_m,Ux_m_s,Uy_m_s,qx_left_m2_s\n";
        forAll(cells, proc) forAll(cells[proc], i)
        {
            out << iteration;
            for (const scalar v : cells[proc][i]) out << ',' << v;
            out << '\n';
        }
        out.close();
        if (!out) FatalErrorInFunction << "Writing flow_cells.csv failed" << exit(FatalError);
    }
}
}
