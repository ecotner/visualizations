#include <math.h>
#include <stdio.h>
#include <string>
#include <stdexcept>

#define RE 0
#define IM 1

enum BC {dirichlet, periodic};

class WaveFunction {
    public:
    int Lx, Ly;

    // Constructor for wavefunction class
    // @param Lx (int): number of grid points in x direction
    // @param Ly (int): number of grid points in y direction
    // @param x0 (double): initial x coordinate of wave packet
    // @param y0 (double): initial y coordinate of wave packet
    // @param s0 (double): initial gaussian wave packet std deviation
    // @param px0 (double): initial momentum in x direction
    // @param py0 (double): initial momentum in y direction
    // @param bdryCond (string): the boundary condition type to use (supports "dirichlet" or "periodic")
    WaveFunction(int Lx, int Ly, double x0, double y0, double s0, double px0, double py0, std::string bdryCond) {
        this->Lx = Lx;
        this->Ly = Ly;
        // allocate memory to psi arrays
        for (int i=0; i<3; i++) {
            psi[i][RE] = (double*) malloc(sizeof(double) * Lx * Ly);
            psi[i][IM] = (double*) malloc(sizeof(double) * Lx * Ly);
        }
        // initialize "current" wave function psi[1]
        if (bdryCond == "dirichlet") {
            this->bdryCond = dirichlet;
            setPsiDirichlet(x0, y0, s0, px0, py0);
        }
        else if (bdryCond == "periodic") {
            this->bdryCond = periodic;
            setPsiPeriodic(x0, y0, s0, px0, py0);
        }
        else throw std::invalid_argument("invalid boundary condition specification: \"" + bdryCond + "\"\n");
    }

    // Destructor; basically just frees memory buffers for psi arrays
    ~WaveFunction() {
        for (int i=0; i<3; i++) {
            free(psi[i][RE]);
            free(psi[i][IM]);
        }
    }

    // Sets/initializes gaussian wave function
    void setPsiDirichlet(double x0, double y0, double s0, double px0, double py0) {
        double rho, theta;
        // set interior points
        for (int x=0; x<Lx; x++) {
            for (int y=0; y<Ly; y++) {
                rho = exp(-(pow(x-x0, 2) + pow(y-y0, 2))/(2*pow(s0, 2)));
                theta = px0*x + py0*y;
                psi[1][RE][Lx*x+y] = rho * cos(theta);
                psi[1][IM][Lx*x+y] = rho * sin(theta);
            }
        }
        // set boundary points to zero
        for (int re=0; re<=1; re++) { // real/imaginary parts
            for (int x=0; x<Lx; x++) { // top/bottom boundaries
                psi[1][re][Lx*x+0] = 0;
                psi[1][re][Lx*x+(Ly-1)] = 0;
            }
            for (int y=0; y<Ly; y++) { // left/right boundaries
                psi[1][re][Lx*0+y] = 0;
                psi[1][re][Lx*(Lx-1)+y] = 0;
            }
        }
    }

    void setPsiPeriodic(double x0, double y0, double s0, double px0, double py0) {
        throw std::invalid_argument("setPsiPeriodic not implemented");
    }

    // Performs FTCS update step
    void stepFTCS(double dt, double dx, double m, double qB) {
        // define a bunch of constants to save us from redoing
        // calculations in the loops
        int xy, xp1, xm1, yp1, ym1;
        double dx2 = dx * dx;
        double qB4dx = (qB*dx/4);
        double q2B24dx2 = qB*qB*dx2/4;
        double dt2mdx2 = dt / (2 * m * dx2);
        // iterate over interior points
        for (int x=1; x<Lx-1; x++) {
            for (int y=1; y<Ly-1; y++) {
                // update some commonly used values/indexes for stencil
                double r2 = x*x + y*y;
                xy = Lx*x+y;
                xp1 = Lx*(x+1)+y;
                xm1 = Lx*(x-1)+y;
                yp1 = Lx*x+(y+1);
                ym1 = Lx*x+(y-1);
                // update real part
                psi[2][RE][xy] = (
                    // kinetic terms
                    dt2mdx2 * (
                        - ( // laplacian (\nabla^2 \psi)
                            psi[1][IM][xp1] + psi[1][IM][xm1] + psi[1][IM][yp1]
                            + psi[1][IM][ym1] - 4*psi[1][IM][xy])
                        + qB4dx * ( // convective derivative (A * \nabla \psi)
                            x*(psi[1][RE][yp1] - psi[1][RE][ym1])
                            -y*(psi[1][IM][xp1] - psi[1][RE][xm1]))
                        + q2B24dx2 * r2 * psi[1][IM][xy] // A^2 * \psi
                    )
                    // + 0*psi[1][IM][xy] // potential (V * \psi)
                );
                // update imaginary part
                psi[2][IM][xy] = (
                    // kinetic terms
                    dt2mdx2 * (
                        + ( // laplacian (\nabla^2 \psi)
                            psi[1][RE][xp1] + psi[1][RE][xm1] + psi[1][RE][yp1]
                            + psi[1][RE][ym1] - 4*psi[1][RE][xy])
                        + qB4dx * ( // convective derivative (A * \nabla \psi)
                            x*(psi[1][IM][yp1] - psi[1][IM][ym1])
                            -y*(psi[1][IM][xp1] - psi[1][IM][xm1]))
                        - q2B24dx2 * r2 * psi[1][RE][xy] // A^2 * \psi
                    )
                    // - dt*0*psi[1][RE][xy] // potential (V * \psi)
                );
            }
        }
        // enforce boundary conditions
        if (bdryCond == dirichlet) {
            // set boundary points to zero
            for (int re=0; re<=1; re++) { // real/imaginary parts
                for (int x=0; x<Lx; x++) { // top/bottom boundaries
                    psi[1][re][Lx*x+0] = 0;
                    psi[1][re][Lx*x+(Ly-1)] = 0;
                }
                for (int y=0; y<Ly; y++) { // left/right boundaries
                    psi[1][re][Lx*0+y] = 0;
                    psi[1][re][Lx*(Lx-1)+y] = 0;
                }
            }
        } else if (bdryCond == periodic) {
            throw std::invalid_argument("FTCS step not implemented for periodic BC");
        }
        // swap array pointers
        for (int i=0; i<=1; i++) {
            double* temp = psi[0][i];
            psi[0][i] = psi[1][i];
            psi[1][i] = psi[2][i];
            psi[2][i] = temp;
        }
    }

    // CTCS update step
    void step() {

    }

    // Return real part of wave function
    double re(double x, double y) {

    }

    // Return imaginary part of wave function
    double im(double x, double y) {

    }

    // Return absolute square of wave function
    double abs2(double x, double y) {

    }

    double abs2(int x, int y) {
        int xy = Lx*x+y;
        return psi[1][RE][xy]*psi[1][RE][xy] + psi[1][IM][xy]*psi[1][IM][xy];
    }

    private:
    // array of pointers to arrays... shape = (# temporal buffers, #re/im components)
    double* psi[3][2];
    int bdryCond;
};

int main() {
    WaveFunction wf = WaveFunction(100, 100, 50.0, 50.0, 5.0, 0.0, 0.0, "dirichlet");
    for (int i=0; i<10; i++) {
        wf.stepFTCS(0.5, 1.0, 1.0, 1e-3);
        printf("%f\n", wf.abs2(50, 50));
    }
}