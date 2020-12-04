#include <math.h>
#include <stdio.h>
#include <string>
#include <stdexcept>

#ifdef __EMSCRIPTEN__
#include <emscripten/bind.h>
#endif

#define RE 0
#define IM 1

enum BC {dirichlet, periodic};

// bilinear interpolation on a square
double interp2d(
    int xL, int xH,
    int yL, int yH,
    double z11, double z12, double z21, double z22,
    double x, double y) {
    // interpolate along x direction first
    double fx1, fx2, dx = xH - xL;
    if (xL == xH) {fx1 = z11; fx2 = z12;}
    else {
        fx1 = z11 + (z21 - z11) * (x - xL) / dx;
        fx2 = z12 + (z22 - z12) * (x - xL) / dx;
    }
    // interpolate along y direction now
    if (yH == yL) return (fx1 + fx2) / 2;
    else return fx1 + (fx2 - fx1) * (y - yL) / (yH - yL);
}

class WaveFunction {
    public:
    int Lx, Ly;
    double dx;

    // Constructor for wavefunction class
    // @param Lx (int): number of grid points in x direction
    // @param Ly (int): number of grid points in y direction
    // @param x0 (double): initial x coordinate of wave packet
    // @param y0 (double): initial y coordinate of wave packet
    // @param s0 (double): initial gaussian wave packet std deviation
    // @param px0 (double): initial momentum in x direction
    // @param py0 (double): initial momentum in y direction
    // @param bdryCond (string): the boundary condition type to use (supports "dirichlet" or "periodic")
    WaveFunction(
        int Lx, int Ly, double dx,
        double x0, double y0,
        double px0, double py0,
        double s0,
        std::string bdryCond) {
        this->Lx = Lx;
        this->Ly = Ly;
        this->dx = dx;
        // allocate memory to psi arrays
        for (int i=0; i<3; i++) {
            psi[i][RE] = (double*) malloc(sizeof(double) * Lx * Ly);
            psi[i][IM] = (double*) malloc(sizeof(double) * Lx * Ly);
        }
        // initialize "current" wave function psi[1]
        if (bdryCond == "dirichlet") {
            this->bdryCond = dirichlet;
            setPsiDirichlet(x0, y0, px0, py0, s0);
        }
        else if (bdryCond == "periodic") {
            this->bdryCond = periodic;
            setPsiPeriodic(x0, y0, px0, py0, s0);
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
    void setPsiDirichlet(double x0, double y0, double px0, double py0, double s0) {
        double rho, theta;
        // set interior points
        for (int x=0; x<Lx; x++) {
            for (int y=0; y<Ly; y++) {
                rho = exp(-(pow(dx*x-x0, 2) + pow(dx*y-y0, 2))/(2*pow(s0, 2)));
                theta = dx*(px0*x + py0*y);
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

    void setPsiPeriodic(double x0, double y0, double px0, double py0, double s0) {
        throw std::invalid_argument("setPsiPeriodic not implemented");
    }

    // Performs FTCS update step
    void stepFTCS(double dt, double B) {
        // define a bunch of constants to save us from redoing
        // calculations in the loops
        int xy, xp1, xm1, yp1, ym1;
        double dx2 = dx * dx;
        double Bdx = B*dx;
        double B2dx4_4 = B*B*dx2*dx2/4;
        double dtdx2 = dt/dx2;
        // iterate over interior points
        for (int x=1; x<Lx-1; x++) {
            for (int y=1; y<Ly-1; y++) {
                // update some commonly used values/indexes for stencil
                xy = Lx*x+y;
                xp1 = Lx*(x+1)+y;
                xm1 = Lx*(x-1)+y;
                yp1 = Lx*x+(y+1);
                ym1 = Lx*x+(y-1);
                // the effective scalar potential
                double Gamma = B2dx4_4 * (x*x + y*y);
                // update real part
                psi[2][RE][xy] = psi[1][RE][xy] + dtdx2 * (
                    - ( // laplacian (\nabla^2 \psi)
                        psi[1][IM][xp1] + psi[1][IM][xm1] + psi[1][IM][yp1]
                        + psi[1][IM][ym1] - 4*psi[1][IM][xy])
                    + Bdx * ( // convective derivative (A * \nabla \psi)
                        x*(psi[1][RE][yp1] - psi[1][RE][ym1])
                        -y*(psi[1][IM][xp1] - psi[1][RE][xm1]))
                    + Gamma * psi[1][IM][xy] // Γ * \psi
                );
                // update imaginary part
                psi[2][IM][xy] = psi[1][IM][xy] + dtdx2 * (
                    + ( // laplacian (\nabla^2 \psi)
                        psi[1][RE][xp1] + psi[1][RE][xm1] + psi[1][RE][yp1]
                        + psi[1][RE][ym1] - 4*psi[1][RE][xy])
                    + Bdx * ( // convective derivative (A * \nabla \psi)
                        x*(psi[1][IM][yp1] - psi[1][IM][ym1])
                        -y*(psi[1][IM][xp1] - psi[1][IM][xm1]))
                    - Gamma * psi[1][RE][xy] // Γ * \psi
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
    void step(double dt, double B) {
        // define a bunch of constants to save us from redoing
        // calculations in the loops
        int xy, xp1, xm1, yp1, ym1;
        double dx2 = dx * dx;
        double Bdx = B*dx;
        double B2dx4_4 = B*B*dx2*dx2/4;
        double dt_dx2 = 2*dt/dx2;
        // iterate over interior points
        for (int x=1; x<Lx-1; x++) {
            for (int y=1; y<Ly-1; y++) {
                // update some commonly used values/indexes for stencil
                double Gamma = B2dx4_4 * (x*x + y*y);
                xy = Lx*x+y;
                xp1 = Lx*(x+1)+y;
                xm1 = Lx*(x-1)+y;
                yp1 = Lx*x+(y+1);
                ym1 = Lx*x+(y-1);
                // update real part
                psi[2][RE][xy] = psi[0][RE][xy] + dt_dx2 * (
                    - ( // laplacian (\nabla^2 \psi)
                        psi[1][IM][xp1] + psi[1][IM][xm1] + psi[1][IM][yp1]
                        + psi[1][IM][ym1] - 4*psi[1][IM][xy])
                    + Bdx * ( // convective derivative (A * \nabla \psi)
                        x*(psi[1][RE][yp1] - psi[1][RE][ym1])
                        -y*(psi[1][IM][xp1] - psi[1][RE][xm1]))
                    + Gamma * psi[1][IM][xy] // A^2 * \psi
                );
                // update imaginary part
                psi[2][IM][xy] = psi[0][IM][xy] + dt_dx2 * (
                    + ( // laplacian (\nabla^2 \psi)
                        psi[1][RE][xp1] + psi[1][RE][xm1] + psi[1][RE][yp1]
                        + psi[1][RE][ym1] - 4*psi[1][RE][xy])
                    + Bdx * ( // convective derivative (A * \nabla \psi)
                        x*(psi[1][IM][yp1] - psi[1][IM][ym1])
                        -y*(psi[1][IM][xp1] - psi[1][IM][xm1]))
                    - Gamma * psi[1][RE][xy] // A^2 * \psi
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

    void multiStep(int n_steps, double dt, double B) {
        for (int i=0; i<n_steps; i++) {
            step(dt, B);
        }
    }

    // Return real part of wave function
    double re(double x, double y) {

    }

    // Return imaginary part of wave function
    double im(double x, double y) {

    }

    // Return absolute square of wave function
    double abs2(double x, double y) {
        int xL = (int) floor(x);
        int xH = (int) ceil(x);
        int yL = (int) floor(y);
        int yH = (int) ceil(y);
        int i11 = Lx*xL+yL;
        int i21 = Lx*xH+yL;
        int i12 = Lx*xL+yH;
        int i22 = Lx*xH+yH;
        double z11 = pow(psi[1][RE][i11], 2) + pow(psi[1][IM][i11], 2);
        double z21 = pow(psi[1][RE][i21], 2) + pow(psi[1][IM][i21], 2);
        double z12 = pow(psi[1][RE][i12], 2) + pow(psi[1][IM][i12], 2);
        double z22 = pow(psi[1][RE][i22], 2) + pow(psi[1][IM][i22], 2);
        return interp2d(xL, xH, yL, yH, z11, z12, z21, z22, x, y);
    }

    // double abs2(int x, int y) {
    //     int xy = Lx*x+y;
    //     return psi[1][RE][xy]*psi[1][RE][xy] + psi[1][IM][xy]*psi[1][IM][xy];
    // }

    int getLx() const {return Lx;}
    int getLy() const {return Ly;}
    double getDx() const {return dx;}

    private:
    // array of pointers to arrays... shape = (# temporal buffers, #re/im components)
    double* psi[3][2];
    int bdryCond;
};

#ifdef __EMSCRIPTEN__
EMSCRIPTEN_BINDINGS(wave_function) {
    emscripten::class_<WaveFunction>("WaveFunction")
    .constructor<int,int,double,double,double,double,double,double,std::string>()
    .function("stepFTCS", &WaveFunction::stepFTCS)
    .function("step", &WaveFunction::step)
    .function("multiStep", &WaveFunction::multiStep)
    .function("abs2", &WaveFunction::abs2)
    .property("Lx", &WaveFunction::getLx)
    .property("Ly", &WaveFunction::getLy)
    .property("dx", &WaveFunction::getDx)
    ;
}

#else
int main() {
    WaveFunction wf = WaveFunction(100, 100, 1.0, 50.0, 50.0, 0.0, 0.0, 5.0, "dirichlet");
    for (int i=0; i<10; i++) {
        wf.stepFTCS(0.5, 1e-3);
        printf("%f\n", wf.abs2(50, 50));
    }
}
#endif