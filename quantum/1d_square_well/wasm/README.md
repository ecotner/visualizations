# Numerical Simulation using WebAssembly

Been learning about the power of WebAssembly (wasm) lately and it seems like the perfect fit for doing relatively simple client-side numerical simulations.
This code here is an example of using wasm to create a simulation of a quantum particle in a box using a finite-difference approximation to the Schrodinger equation for an infinite square well potential.
The mechanics of the simulation are written in C++ (`wavefunction.cpp` file), where we define a `WaveFunction` class which initializes the wavefunction in the class constructor and evolves the state through time using a `.step()` method.
We can use the Emscripten compiler `emcc` to produce a `module.wasm` WebAssembly binary and `module.js` "glue code" in order to access and manipulate the `WaveFunction` object directly in JavaScript (the `script.js` file)!

[Click here](./index.html) to see the simulation in action!