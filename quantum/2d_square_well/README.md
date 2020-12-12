# 2D quantum simulation in electromagnetic potential

Simulation of a particle in a 2D box in the presence of electromagnetic fields. [See the simulation in action!](https://ecotner.github.io/visualizations/quantum/2d_square_well/)

## Physics
This simulation basically solves the Schrodinger equation coupled to a uniform electric/magnetic field configuration.
It utilizes an [interesting algorithm](https://aip.scitation.org/doi/pdf/10.1063/1.168483) that is both an _explicit_ numerical scheme and unconditionally stable thanks to a clever decomposition of the Hamiltonian into groups of commuting two-state systems whose time evolution can be solved for analytically.
This allows for probability mass to be explicitly conserved (within rounding error), so that you can take time steps of arbitrary magnitude without getting exploding solutions.
I have outlined some of the mathematical details of the algorithm (including some corrections of errors in the original paper!) in [this jupyter notebook](math.ipynb) if you would like to read more.
(You can also enjoy my `mspaint` diagrams trying to explain commutivity of creation/annihilation operators on a lattice!)

## Technology
This simulation is able to run inside your browser (or even on your phone) thanks to the magic of _WebAssembly_ (aka WASM); a web standard for compiled binaries that is very performant, and can be compiled from many different languages, including C, C++, Rust, Java, and JavaScript.
This specific project utilizes a C/C++ -> WASM compiler called [`emscripten`](https://emscripten.org/index.html); particularly a library called [`embind`](https://emscripten.org/docs/porting/connecting_cpp_and_javascript/embind.html) that is able to create JavaScript bindings for C++ classes.
A `WaveFunction` class is defined in the `wavefunction_unitary.cpp` file which contains the complete state of the quantum wavefunction and exposes methods for evolving it in time, measuring the probability density/mass, etc.

While the numerical heavy-lifting is taken care of by WASM, the JavaScript acts as a "glue" that links the front-end to the wave function object.
Numerical results are accessed on the JS side and rendered using a canvas HTML element after each time step, and sliders allow the user to tweak parameters and run/pause/restart the simulation on demand.

## Compliation/running
If you want to get this working on your local machine, you can fisrt compile the WASM binary from the C++ source by installing the `emscripten` SDK and running
```
emcc --bind -o module.js wavefunction_unitary.cpp -O3 -s WASM=1 -s ALLOW_MEMORY_GROWTH=1
```
Then you can start a local webserver (WASM will not run as a static file)
```
python3 -m http.server 8000
```

## TODO
* Split the wave function into real/imaginary components and get rid of the `complex<double>` types to improve performance.
* Find gauge-invariant way to specify initial wave packet momentum? Might not be possible!