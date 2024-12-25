#include "inmost.h"
#include <stdio.h>
#include <math.h>

using namespace INMOST;
using namespace std;

const double dx = 1.0;
const double dy = 1.0;
const double dxy = 0.0;
const double a = 10;

double C(double x, double y) // analytical solution
{
    return sin(a*x) * sin(a*y);
}

double source(double x, double y)
{
    return -a*a * (2.*dxy * cos(a*x)*cos(a*y) - (dx+dy) * sin(a*x)*sin(a*y));
}

enum BoundCondType
{
    BC_DIR = 1,
    BC_NEUM = 2
};

// Class including everything needed
class Problem
{
private:
    /// Mesh
    Mesh &m;
    // =========== Tags =============
    /// Solution tag: 1 real value per cell
    Tag tagConc;
    /// Diffusion tensor tag: 3 real values (Dx, Dy, Dxy) per cell
    Tag tagD;
    /// Boundary condition type tag: 1 integer value per face, sparse on faces
    Tag tagBCtype;
    /// Boundary condition value tag: 1 real value per face, sparse on faces
    Tag tagBCval;
    /// Right-hand side tag: 1 real value per cell
    Tag tagSource;
    /// Analytical solution tag: 1 real value per cell
    Tag tagConcAn;
    /// Global index tag: 1 integer value per cell
    Tag tagGlobInd;
    /// Boundary conductivity: 1 real value per face, sparse on faces
    Tag tagBCcond;

    // =========== Tag names ===========
    const string tagNameConc = "Concentration";
    const string tagNameD = "Diffusion_tensor";
    const string tagNameBCtype = "BC_type";
    const string tagNameBCval = "BC_value";
    const string tagNameSource = "Source";
    const string tagNameConcAn = "Concentration_analytical";
    const string tagNameGlobInd = "Global_Index";
    const string tagNameBCcond = "BC_conductivity";

public:
    Problem(Mesh &m_);
    ~Problem();
    void initProblem();
    void assembleGlobalSystem(Sparse::Matrix &A, Sparse::Vector &rhs);
    void run();
};

Problem::Problem(Mesh &m_) : m(m_)
{
}

Problem::~Problem()
{
}


double calc_tf(rMatrix const& D, rMatrix const& nf, double *dA){
    rMatrix DdA(2, 1);
    DdA(0, 0) = D(0, 0) * dA[0] + D(0, 1) * dA[1];
    DdA(1, 0) = D(1, 0) * dA[0] + D(1, 1) * dA[1];
    double temp = (DdA(0, 0) * nf(0, 0) + DdA(1, 0) * nf(1, 0)) /
                    (dA[0]* dA[0] + dA[1] * dA[1]);
    return temp;
}


void Problem::initProblem()
{
    // Init tags
    tagConc = m.CreateTag(tagNameConc, DATA_REAL, CELL, NONE, 1);
    tagD = m.CreateTag(tagNameD, DATA_REAL, CELL, NONE, 3);
    tagBCtype = m.CreateTag(tagNameBCtype, DATA_INTEGER, FACE, FACE, 1);
    tagBCval = m.CreateTag(tagNameBCval, DATA_REAL, FACE, FACE, 1);
    tagSource = m.CreateTag(tagNameSource, DATA_REAL, CELL, CELL, 1);
    tagConcAn = m.CreateTag(tagNameConcAn, DATA_REAL, CELL, CELL, 1);
    tagGlobInd = m.CreateTag(tagNameGlobInd, DATA_INTEGER, CELL, NONE, 1);
    tagBCcond = m.CreateTag(tagNameBCcond, DATA_REAL, FACE, FACE, 1);

    // Cell loop
    // 1. Set diffusion tensor values
    // 2. Write analytical solution and source tags
    // 3. Assign global indices
    int glob_ind = 0;
    for(Mesh::iteratorCell icell = m.BeginCell(); icell != m.EndCell(); icell++){
        Cell c = icell->getAsCell();
        c.RealArray(tagD)[0] = dx; // Dx
        c.RealArray(tagD)[1] = dy; // Dy
        c.RealArray(tagD)[2] = dxy; // Dxy
        double xc[2];
        c.Barycenter(xc);
        c.Real(tagConcAn) = C(xc[0], xc[1]);
        c.Real(tagSource) = source(xc[0], xc[1]);
        c.Integer(tagGlobInd) = glob_ind;
        glob_ind++;
    }

    // Face loop:
    // 1. Set BC
    for(Mesh::iteratorFace iface = m.BeginFace(); iface != m.EndFace(); iface++){
        Face f = iface->getAsFace();
        double xf[2];
        f.Barycenter(xf);
        if(f.Boundary()) {
            f.Integer(tagBCtype) = BC_DIR;
            f.Real(tagBCval) = C(xf[0], xf[1]);
        } else {
            Cell cA, cB;
            cA = f.BackCell();
            cB = f.FrontCell();
            if(!cB.isValid()){
                std::cout << "Invalid FrontCell!" << endl;
                exit(1);
            }

            double xA[2], xB[2];
            cA.Barycenter(xA), cB.Barycenter(xB);
            rMatrix nf(2,1);
            f.UnitNormal(nf.data());
            double dA[2], dB[2];
            for (int i = 0; i < 2; i++) {
                dA[i] = xf[i] - xA[i];
                dB[i] = xf[i] - xB[i];
            }

            rMatrix DA(2, 2);
            DA(0, 0) = cA.RealArray(tagD)[0];
            DA(0, 1) = cA.RealArray(tagD)[2];
            DA(1, 0) = cA.RealArray(tagD)[2];
            DA(1, 1) = cA.RealArray(tagD)[1];
            double tfA = calc_tf(DA, nf, dA);

            rMatrix DB(2, 2);
            DB(0, 0) = cB.RealArray(tagD)[0];
            DB(0, 1) = cB.RealArray(tagD)[2];
            DB(1, 0) = cB.RealArray(tagD)[2];
            DB(1, 1) = cB.RealArray(tagD)[1];
            double tfB = calc_tf(DB, nf, dB);

            f.Real(tagBCcond) = tfA * tfB / (tfA - tfB);
        }
    }
}

void Problem::assembleGlobalSystem(Sparse::Matrix &M, Sparse::Vector &rhs)
{
    // Face loop
    // Calculate transmissibilities using
    // two-point flux approximation (TPFA)
    for(Mesh::iteratorFace iface = m.BeginFace(); iface != m.EndFace(); iface++){
        Face f = iface->getAsFace();
        double xf[2];
        rMatrix nf(2,1);
        f.UnitNormal(nf.data());
        f.Barycenter(xf);
        if(f.Boundary()){
            int BCtype = f.Integer(tagBCtype);
            if(BCtype == BC_NEUM){
                continue;
            }
            else if(BCtype == BC_DIR){
                // implement by yourself
                Cell A;
                A = f.BackCell();

                double xA[2];
                A.Barycenter(xA);

                double dA[2] = {xf[0] - xA[0], xf[1] - xA[1]};

                rMatrix DA(2, 2);
                DA(0, 0) = A.RealArray(tagD)[0];
                DA(0, 1) = A.RealArray(tagD)[2];
                DA(1, 0) = A.RealArray(tagD)[2];
                DA(1, 1) = A.RealArray(tagD)[1];

                double t = calc_tf(DA, nf, dA); // transmissibility
                
                int id = A.Integer(tagGlobInd);
                M[id][id] -= t * f.Area();
                rhs[id] -= t * f.Real(tagBCval) * f.Area();
            }
        }
        else{
            // Internal face
            Cell cA, cB;
            cA = f.BackCell();
            cB = f.FrontCell();
            if(!cB.isValid()){
                cout << "Invalid FrontCell!" << endl;
                exit(1);
            }

            // implement by yourself
            double t = f.Real(tagBCcond);
            int idA = cA.Integer(tagGlobInd);
            int idB = cB.Integer(tagGlobInd);
            M[idA][idA] += t * f.Area();
            M[idA][idB] -= t * f.Area();
            M[idB][idA] -= t * f.Area();
            M[idB][idB] += t * f.Area();
        }
    }
    for(auto icell = m.BeginCell(); icell != m.EndCell(); icell++){
        Cell c = icell->getAsCell();
        int i = c.Integer(tagGlobInd);
        rhs[i] -= c.Real(tagSource) * c.Volume();
    }
}

void Problem::run()
{
    // Matrix size
    unsigned N = static_cast<unsigned>(m.NumberOfCells());
    // Global matrix called 'stiffness matrix'
    Sparse::Matrix A;
    // Solution vector
    Sparse::Vector sol;
    // Right-hand side vector
    Sparse::Vector rhs;
    std::cout << "N = " << N << "\n";

    A.SetInterval(0, N);
    sol.SetInterval(0, N);
    rhs.SetInterval(0, N);

    assembleGlobalSystem(A, rhs);

    string solver_name = "inner_mptiluc";
    Solver S(solver_name);
    S.SetParameter("drop_tolerance", "0");
    S.SetParameter("absolute_tolerance", "1e-14");
    S.SetParameter("relative_tolerance", "1e-10");


    S.SetMatrix(A);
    bool solved = S.Solve(rhs, sol);
    printf("Number of iterations: %d\n", S.Iterations());
    printf("Residual:             %e\n", S.Residual());
    if(!solved){
        printf("Linear solver failed: %s\n", S.GetReason().c_str());
        exit(1);
    }

    double normC = 0.0, normL2 = 0.0;
    for(Mesh::iteratorCell icell = m.BeginCell(); icell != m.EndCell(); icell++){
        Cell c = icell->getAsCell();
        unsigned ind = static_cast<unsigned>(c.Integer(tagGlobInd));
        c.Real(tagConc) = sol[ind];
        double diff = fabs(c.Real(tagConc) - c.Real(tagConcAn));
        normL2 += diff * c.Volume();
        normC = max(normC, diff);
    }
    printf("\nError C-norm:  %e\n", normC);
    printf("Error L2-norm: %e\n", normL2);

    m.Save("res.pvtk");
}

int main(int argc, char ** argv)
{
    if( argc < 2 )
    {
        printf("Usage: %s mesh_file\n", argv[0]);
        return -1;
    }
    for (int i = 1; i < argc; i++) {
        Mesh m;
        m.Load(argv[i]);
        Problem P(m);
        P.initProblem();
        P.run();
        printf("Success\n\n");
    }
    return 0;
}
