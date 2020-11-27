// Class for defining wavefunction and time evolution of a quantum
// particle in a 1D infinite square well

#ifdef __EMSCRIPTEN__
#include <emscripten/bind.h>
#endif
#include <stdlib.h>
#include <string>
#include <math.h>
#include <algorithm>


class WaveFunction {
    public:
    double *psi_next; // just a buffer; do not access directly
    double *psi;
    double *psi_prev;

    // Constructor; initializes wave function values
    // @param num_points (int): number of grid points to simulate
    // @param p0 (double): initial momentum, as multiples of 2*pi/(L+1)
    // @returns NULL
    WaveFunction(int length, double dt_, double p0) {
        L = length;
        dt = dt_;
        // psi will contain both real/imag parts; alternate elements
        // starting with psi[0] = real, psi[1] = imag, ...
        // Will initialize three arrays to hold everything
        psi_next = (double*) malloc(sizeof(double) * 2 * (L+1));
        psi = (double*) malloc(sizeof(double) * 2 * (L+1));
        psi_prev = (double*) malloc(sizeof(double) * 2 * (L+1));
        resetPsi(2.0*M_PI*p0/L, "gaussian");
        // halfStep();
        for (int i=0; i<L+1; i++) {
            psi_prev[2*i] = psi[2*i];
            psi_prev[2*i+1] = psi[2*i+1];
        }
        dt /= 1000;
        multiStep(1000);
        dt *= 1000;
    }

    // destructor; frees all the memory contained in the psi arrays
    ~WaveFunction() {
        free(psi_next);
        free(psi);
        free(psi_prev);
    }

    // resets/initializes the wave function
    // @param p0 (double): the initial momentum
    // @param type (string): the type of initial wave function; right now "gaussian" is the only option
    void resetPsi(double p0, std::string type) {
        double rho;
        double sigma = L/10.;
        double p = 0.25;
        // double N;
        for (int x=0; x<L+1; x++) {
            if (type == "gaussian") {
                rho = exp(-pow(x - L/2, 2)/(2*pow(sigma, 2)));
            }
            // rho *= (pow(x, p) * pow(L-1-x, p)) / pow((L-1)/2, 2*p);
            psi[2*x] = rho * cos(p0*x);
            psi[2*x+1] = rho * sin(p0*x);
        }
        psi[0] = 0.; psi[1] = 0.;
        psi[2*L] = 0.; psi[2*L+1] = 0.;
    }

    // Performs a FTCS update; this is numerically unstable long-term, but
    // will work to initialize the first two timesteps
    void halfStep() {
        // iterate over all non-boundary points
        for (int x=1; x<L; x++) {
            psi_next[2*x] = psi[2*x] + 0.5 * dt * (psi[2*(x-1)+1] - 2*psi[2*x+1] + 2*psi[2*(x+1)+1]);
            psi_next[2*x+1] = psi[2*x+1] - 0.5 * dt * (psi[2*(x-1)] - 2*psi[2*x] + 2*psi[2*(x+1)]);
        }
        // update boundary points
        psi_next[0] = 0.; psi_next[1] = 0.; // left side
        psi_next[2*L] = 0.; psi_next[2*L+1] = 0.; // right side
        // psi_next[0] = psi[0]; psi_next[1] = psi[1]; // left side
        // psi_next[2*L] = psi[2*L]; psi_next[2*L+1] = psi[2*L+1]; // right side
        // swap pointers
        double* temp = psi_prev;
        this->psi_prev = psi;
        this->psi = psi_next;
        this->psi_next = temp;
    }

    // Performs a single CTCS update; should be stable numerically for dt < 1
    void step() {
        // iterate over all non-boundary points
        double lap_re, lap_im;
        for (int x=1; x<L; x++) {
            lap_re = psi[2*(x-1)] - 2*psi[2*x] + 2*psi[2*(x+1)];
            lap_im = psi[2*(x-1)+1] - 2*psi[2*x+1] + 2*psi[2*(x+1)+1]; // creates problems on right boundary
            psi_next[2*x] = psi_prev[2*x] + dt * lap_im;
            psi_next[2*x+1] = psi_prev[2*x+1] - dt * lap_re;
        }
        // update boundary points
        psi_next[0] = 0.; psi_next[1] = 0.; // left side
        psi_next[2*L] = 0.; psi_next[2*L+1] = 0.; // right side
        // psi_next[0] = psi[0]; psi_next[1] = psi[1]; // left side
        // psi_next[2*L] = psi[2*L]; psi_next[2*L+1] = psi[2*L+1]; // right side
        // swap pointers
        double* temp = psi_prev;
        this->psi_prev = psi;
        this->psi = psi_next;
        this->psi_next = temp;
    }

    // Convenience to just do multiple steps at once
    // @param num_steps (int): number of update steps to do at once
    void multiStep(int num_steps) {
        for (int i=0; i<num_steps; i++) {
            step();
        }
    }

    // access (interpolated) wave function values

    // get real part of wave function
    double getPsiRe(double x) {
        // get two high/low bounding grid points
        int xH = (int) ceil(x);
        int xL = (int) floor(x);
        // interpolate between the two points
        if (xH == xL) {
            return psi[2*xL];
        }
        double yH = psi[2*xH], yL = psi[2*xL];
        return yL + ((yH-yL)/(xH-xL)) * (x-xL);
    }

    // get imaginary part of wave function
    double getPsiIm(double x) {
        // get two high/low bounding grid points
        int xH = (int) ceil(x);
        int xL = (int) floor(x);
        // interpolate between the two points
        if (xH == xL) {
            return psi[2*xL+1];
        }
        double yH = psi[2*xH+1], yL = psi[2*xL+1];
        return yL + ((yH-yL)/(xH-xL)) * (x-xL);
    }

    // get absolute square of wave function
    double getPsiAbsSq(double x) {
        // get two high/low bounding grid points
        int xH = (int) ceil(x);
        int xL = (int) floor(x);
        // interpolate between the two points
        if (xH == xL) {
            return pow(psi[2*xL], 2) + pow(psi[2*xL+1], 2);
        }
        double yH = pow(psi[2*xH], 2) + pow(psi[2*xH+1], 2);
        double yL = pow(psi[2*xL], 2) + pow(psi[2*xL+1], 2);
        return yL + ((yH-yL)/(xH-xL)) * (x-xL);
    }

    int getL() {
        return L;
    }

    private:
    int L;
    double dt;
};


int main() {
    WaveFunction psi(500, 1e-3, 2);
    // printf("wave function values:\n");
    // int L = psi.getL();
    // for (int i=0; i<L; i++) {
    //     printf("(%f, %f)^2 = %f\n", psi.getPsiRe(i), psi.getPsiIm(i), psi.getPsiAbsSq(i));
    // }
    for (int i=0; i<4; i++) {
        printf("prev: %p\n", psi.psi_prev);
        printf("cur: %p\n", psi.psi);
        printf("next: %p\n\n", psi.psi_next);
        psi.step();
    }
}

#ifdef __EMSCRIPTEN__
EMSCRIPTEN_BINDINGS(wave_function_module) {
    emscripten::class_<WaveFunction>("WaveFunction")
    .constructor<int, double, double>()
    .function("halfStep", &WaveFunction::halfStep)
    .function("step", &WaveFunction::step)
    .function("multiStep", &WaveFunction::multiStep)
    .function("getPsiRe", &WaveFunction::getPsiRe)
    .function("getPsiIm", &WaveFunction::getPsiIm)
    .function("getPsiAbsSq", &WaveFunction::getPsiAbsSq)
    .function("getL", &WaveFunction::getL)
    ;
}
#endif