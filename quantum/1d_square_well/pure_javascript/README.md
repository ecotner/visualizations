# Client-side numerical computations utilizing WebGL

This simple script utilizes the `gpu.js` project to perform _client-side GPGPU_ calculations in order to simulate the 1D Schrodinger equation.
`gpu.js` allows you to write GPU kernels using pure JavaScript, which is then transpiled into WebGL.
Some care needs to be taken when writing efficient GPU kernels to avoid warp divergences due to conditional branching, so if you look at the code and ask "why is he making floating-point versions of inequality operators?", that's why.

[Check out](https://ecotner.github.io/visualizations/quantum/1d_square_well/pure_javascript/) the simulation in action.
You might need to go into your browser's settings and enable WebGL or "hardware acceleration" to be able to see it.
