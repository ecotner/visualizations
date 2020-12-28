# Fractal explorer

Simple javascript implementation of an interactive fractal exploration tool.
When you scroll over the canvas, a new frame is computed which iterates over each pixel and calculates the number of iterations of the update step z <- z^2 + c needed until |z| > 2.
Some unique color is then assigned to each integer and is used to draw the frame.
[Click here](https://ecotner.github.io/visualizations/fractals/single-thread/) to see the explorer in action!
