# Fractal explorers

This project is a set of live fractal renderers that calculate the number of escape iterations under the map z <- z^2 + c.
I have two different implementations:

* One is a single-threaded CPU application written in pure javascript. This one can get down to about 10^13 zoom due to use of 64-bit precision floating point.
* The other one uses `gpu.js` (which utilizes WebGL) to parallelize over each pixel being rendered. It is _much_ faster than the other one, but can only get down to about 10^6 zoom due to only being able to use 32-bit floating point precision