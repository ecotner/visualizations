# 2D quantum simulation in electromagnetic potential

A particle in a 2D box in the presence of electromagnetic fields (just magnetic for now).

Found an interesting algorithm based on analytic calculation of of the time evolution operator (with an approximate hamiltonian) [here](https://aip.scitation.org/doi/pdf/10.1063/1.168483).

To compile from source, run
```
emcc --bind -o module.js wavefunction_unitary.cpp -Os -s WASM=1
```

### TODO:
* it seems like the y direction is swapped for some reason
    * when specifying direction of initial momentum, angle appears to work in clockwise direction
    * when specifying direction of electric field, seems to be all messed up
    * when specifying initial position of wavefunction, everything appears to work as expected
    * how to rectify all these things???
* wavefunction spreads more in the y direction than it does in the x direction, even in the absence of external fields or with starting momenta