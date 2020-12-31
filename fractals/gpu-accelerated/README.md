# Fractal explorer

Simple javascript implementation of an interactive fractal exploration tool which utilizes [gpu.js](https://github.com/gpujs/gpu.js/) to accelerate the computation of each visualized frame.
When you scroll over the canvas, a new frame is computed which calculates the number of iterations of the update step z <- z^2 + c needed until |z| > 2 for each pixel in parallel by putting them on the GPU and running the iteration loop.
Some unique color is then assigned to each integer and is used to draw the frame.
[Click here](https://ecotner.github.io/visualizations/fractals/gpu-accelerated/index.html) to see the explorer in action!
