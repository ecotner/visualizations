// Author: Eric Cotner
// Date: 2020-12-05
//
// This file's purpose is to create an object-oriented interface for manipulating
// and propagating the wave function for a 2D quantum particle in an electromagnetic
// field using a unitary time evolution operator. The algorithm is based on the
// following paper:
//
// Title: Algorithm to solve the time-dependent Schrödinger equation for a charged
// particle in an inhomogeneous magnetic field: Application to the Aharonov–Bohm effect
// Ref: Computers in Physics 8, 600 (1994); doi: 10.1063/1.168483
// Link: https://aip.scitation.org/doi/pdf/10.1063/1.168483

#include <stdlib.h>
#include <complex>
#include <math.h>
#include <stdio.h>
#ifdef __EMSCRIPTEN__
#include <emscripten/bind.h>
#endif

using namespace std;
complex<double> I (0, 1);

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

    // Wavefunction simulation constructor
    // @param Nx (int): Number of spatial grid points in x direction
    // @param Ny (int): Number of spatial grid points in y direction
    // @param dx (double): Spatial separation between grid points in units of 1/m
    // @param dt (double): Temporal separation between time steps in units of 2/m
    // @param B (double): Magnetic field strength in units of m^2/q
    WaveFunction(
        int Nx,
        int Ny,
        double dx,
        double dt,
        double B,
        double Ex,
        double Ey,
        double x0,
        double y0,
        double s0,
        double px0,
        double py0
    ) {
        // add attributes
        this->Nx = Nx;
        this->Ny = Ny;
        this->dx = dx;
        this->dt = dt;
        this->B = B;
        // Allocate memory for all the dynamically-sized arrays
        psi = (complex<double>*) malloc(sizeof(complex<double>) * Nx * Ny);
        for (int i=0; i<3; i++) {
            ca[i] = (double*) malloc(sizeof(double) * Nx * Ny);
            for (int j=0; j<2; j++) {
                sa[j][i] = (complex<double>*) malloc(sizeof(complex<double>) * Nx * Ny);
            }
        }
        for (int i=0; i<2; i++) {
            cb[i] = (double*) malloc(sizeof(double) * Nx * Ny);
            for (int j=0; j<2; j++) {
                sb[j][i] = (complex<double>*) malloc(sizeof(complex<double>) * Nx * Ny);
            }
        }
        // Precompute/cache repeatedly used variables
        precomputeArrays(B, Ex, Ey);
        // Initialize wave function
        initializePsi(x0, y0, s0, px0, py0);
    }

    // Pre-computes the ca, cb, sa, sb arrays to save FLOPS during time evolution
    void precomputeArrays(double B, double Ex, double Ey) {
        // Assuming uniform magnetic field in z direction with gauge \vec{A} = (A, 0, 0) implies A = -B*y
        double Ca0 = 5./(dx*dx);
        double Ca1 = -(4./(3*dx*dx));
        double Ca2 = 1./(12*dx*dx);
        complex<double> a[3]; // \alpha coefficients of Hamiltonian
        double b[2]; // \beta coefficients of Hamiltonian
        double A0, A1, A2; // Gauge fields
        for (int i=0; i<Nx; i++) {
            for (int j=0; j<Ny; j++) {
                int xy = Nx*i+j;
                A0 = -B*(j*dx); A1 = -B*((j+1)*dx); A2 = -B*((j+2)*dx);
                // alpha computations
                a[0] = Ca0 + A0*A0 - dx*(i*Ex + j*Ey);
                a[1] = Ca1 * (1. - I*0.5*dx*(A0 + A1));
                a[2] = Ca2 * (1. - I*dx*(A0 + A2));
                for (int k=0; k<3; k++) {
                    ca[k][xy] = cos(abs(a[k])*dt);
                    sa[0][k][xy] = I * (a[k]/abs(a[k])) * sin(abs(a[k])*dt);
                    sa[1][k][xy] = I * (conj(a[k])/abs(a[k])) * sin(abs(a[k])*dt);
                }

                // beta computations
                b[0] = 16. * Ca2; b[1] = Ca2;
                for (int k=0; k<2; k++) {
                    cb[k][xy] = cos(abs(b[k])*dt);
                    sb[0][k][xy] = I * (b[k]/abs(b[k])) * sin(abs(b[k])*dt);
                    sb[1][k][xy] = I * (conj(b[k])/abs(b[k])) * sin(abs(b[k])*dt);
                }
            }
        }
    }

    // Initialize the wave function
    void initializePsi(double x0, double y0, double s0, double px0, double py0) {
        double rho, theta;
        // set interior points
        for (int i=0; i<Nx; i++) {
            for (int j=0; j<Ny; j++) {
                rho = exp(-(pow(dx*i-x0, 2) + pow(dx*j-y0, 2))/(2*pow(s0, 2)));
                theta = dx*(px0*i + py0*j);
                psi[Nx*i+j] = rho * complex<double>(cos(theta), sin(theta));
            }
        }
        // Dirichlet boundary conditions are automatically enforced?
    }

    // Time evolution operators for sub-Hamiltonians H1 - H9

    // two-hop; alpha; 4i, 4i+1
    void substep_H1() {
        for (int i=0; i<Nx-3; i+=4) {
            for (int j=0; j<Ny; j++) {
                for (int k=0; k<=1; k++) {
                    int xy = Nx*(i+k)+j;
                    int xp2 = Nx*(i+2+k)+j;
                    complex<double> _psi (psi[xy]);
                    psi[xy] = ca[2][xy] * _psi - sa[0][2][xy] * psi[xp2];
                    psi[xp2] = ca[2][xy] * psi[xp2] - sa[1][2][xy] * _psi;
                }
            }
        }
    }

    // two-hop; alpha; 4i+2, 4i+3
    void substep_H2() {
        for (int i=0; i<Nx-5; i+=4) {
            for (int j=0; j<Ny; j++) {
                for (int k=2; k<=3; k++) {
                    int xy = Nx*(i+k)+j;
                    int xp2 = Nx*(i+2+k)+j;
                    complex<double> _psi (psi[xy]);
                    psi[xy] = ca[2][xy] * _psi - sa[0][2][xy] * psi[xp2];
                    psi[xp2] = ca[2][xy] * psi[xp2] - sa[1][2][xy] * _psi;
                }
            }
        }
    }

    // one-hop; alpha; 2i
    void substep_H3() {
        for (int i=0; i<Nx-1; i+=2) {
            for (int j=0; j<Ny; j++) {
                int xy = Nx*i+j;
                int xp1 = Nx*(i+1)+j;
                complex<double> _psi (psi[xy]);
                psi[xy] = ca[1][xy] * _psi - sa[0][1][xy] * psi[xp1];
                psi[xp1] = ca[1][xy] * psi[xp1] - sa[1][1][xy] * _psi;
            }
        }
    }

    // one-hop; alpha; 2i+1
    void substep_H4() {
        for (int i=0; i<Nx-2; i+=2) {
            for (int j=0; j<Ny; j++) {
                int xy = Nx*(i+1)+j;
                int xp1 = Nx*(i+2)+j;
                complex<double> _psi (psi[xy]);
                psi[xy] = ca[1][xy] * _psi - sa[0][1][xy] * psi[xp1];
                psi[xp1] = ca[1][xy] * psi[xp1] - sa[1][1][xy] * _psi;
            }
        }
    }

    // two-hop; beta; 4j, 4j+1
    void substep_H5() {
        for (int i=0; i<Nx; i++) {
            for (int j=0; j<Ny-3; j+=4) {
                for (int k=0; k<=1; k++) {
                    int xy = Nx*i+(j+k);
                    int yp2 = Nx*i+(j+2+k);
                    complex<double> _psi (psi[xy]);
                    psi[xy] = cb[1][xy] * _psi - sb[0][1][xy] * psi[yp2];
                    psi[yp2] = cb[1][xy] * psi[yp2] - sb[1][1][xy] * _psi;
                }
            }
        }
    }

    // two-hop; beta; 4j+2, 4i+3
    void substep_H6() {
        for (int i=0; i<Nx; i++) {
            for (int j=0; j<Ny-5; j+=4) {
                for (int k=2; k<=3; k++) {
                    int xy = Nx*i+(j+k);
                    int yp2 = Nx*i+(j+2+k);
                    complex<double> _psi (psi[xy]);
                    psi[xy] = cb[1][xy] * _psi - sb[0][1][xy] * psi[yp2];
                    psi[yp2] = cb[1][xy] * psi[yp2] - sb[1][1][xy] * _psi;
                }
            }
        }
    }

    // one-hop; beta; 2j
    void substep_H7() {
        for (int i=0; i<Nx; i++) {
            for (int j=0; j<Ny-1; j+=2) {
                int xy = Nx*i+j;
                int yp1 = Nx*i+(j+1);
                complex<double> _psi (psi[xy]);
                psi[xy] = cb[0][xy] * _psi - sb[0][0][xy] * psi[yp1];
                psi[yp1] = cb[0][xy] * psi[yp1] - sb[1][0][xy] * _psi;
            }
        }
    }

    // one-hop; beta; 2j+1
    void substep_H8() {
        for (int i=0; i<Nx; i++) {
            for (int j=0; j<Ny-2; j+=2) {
                int xy = Nx*i+(j+1);
                int yp1 = Nx*i+(j+2);
                complex<double> _psi (psi[xy]);
                psi[xy] = cb[0][xy] * _psi - sb[0][0][xy] * psi[yp1];
                psi[yp1] = cb[0][xy] * psi[yp1] - sb[1][0][xy] * _psi;
            }
        }
    }

    // zero-hop; alpha
    void substep_H9() {
        for (int i=0; i<Nx; i++) {
            for (int j=0; j<Ny; j++) {
                int xy = Nx*i+j;
                // simple phase rotation
                psi[xy] = (ca[0][xy] - sa[0][0][xy]) * psi[xy];
            }
        }
    }

    // Performs a single time step
    void step(int n=1) {
        for (int i=0; i<n; i++) {
            substep_H1();
            substep_H2();
            substep_H3();
            substep_H4();
            substep_H5();
            substep_H6();
            substep_H7();
            substep_H8();
            substep_H9();
        }
    }

    // Return absolute square of wave function at interpolated point
    double abs2(double x, double y) {
        x /= dx; y /= dx;
        int xL = (int) floor(x);
        int xH = (int) ceil(x);
        int yL = (int) floor(y);
        int yH = (int) ceil(y);
        double z11 = norm(psi[Nx*xL+yL]);
        double z21 = norm(psi[Nx*xH+yL]);
        double z12 = norm(psi[Nx*xL+yH]);
        double z22 = norm(psi[Nx*xH+yH]);
        return interp2d(xL, xH, yL, yH, z11, z12, z21, z22, x, y);
    }

    // Return total probability mass of wave function
    double probMass() {
        double prob = 0;
        for (int i=0; i<Nx; i++) {
            for (int j=0; j<Ny; j++) {
                prob += norm(psi[Nx*i+j]);
            }
        }
        return dx*dx*prob;
    }

    // Destructor; frees memory
    ~WaveFunction() {
        free(psi);
        for (int i=0; i<3; i++) {
            free(ca[i]);
            for (int j=0; j<2; j++) {
                free(sa[j][i]);
            }
        }
        for (int i=0; i<2; i++) {
            free(cb[i]);
            for (int j=0; j<2; j++) {
                free(sb[j][i]);
            }
        }
    }

    private:
    int Nx, Ny; // # of grid points
    double dx, dt, B; // parameters
    complex<double> *psi; // wave function values
    // diagonal cos(|a|t) and cos(|b|t) rotation matrix elements; index is hop size
    double *ca[3], *cb[2];
    // off-diagonal sin(|a|t) and sin(|b|t) rotation matrix elements; first index is
    // row #, second index is hop size; sb=0 for n=0, so there are only two (0, 1) -> (1, 2)
    complex<double> *sa[2][3], *sb[2][2];
};

#ifdef __EMSCRIPTEN__
EMSCRIPTEN_BINDINGS(wave_function) {
    emscripten::class_<WaveFunction>("WaveFunction")
    .constructor<int,int,double,double,double,double,double,double,double,double,double,double>()
    .function("step", &WaveFunction::step)
    .function("abs2", &WaveFunction::abs2)
    .function("probMass", &WaveFunction::probMass)
    ;
}

#else
int main() {
    double dx = 1e-3;
    double dt = 1e-1;
    int N = 100;
    WaveFunction wf (
        N, N, // Nx, Ny
        dx, dt, // dx, dt
        1., 0., 0., // B, Ex, Ey
        5*dx, 5*dx, 2.*dx, // x0, y0, s0
        0.5*2*M_PI/(N*dx), 0.5*2*M_PI/(N*dx) // px0, py0
    );
    printf("prob. mass before step: %f\n", wf.probMass());
    for (int i=0; i<1000; i++) wf.step();
    printf("prob. mass after step: %f\n", wf.probMass());
}

#endif