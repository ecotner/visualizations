// pretty basic live fractal renderer. iterates through every pixel in the canvas
// and calculates $z <- z^2 + c$ until $|z| > R$. then the number of iterations to
// "escape" are given a nice color and plotted. each pixel is calculated sequentially,
// so it's pretty slow. we could speed it up by either using web workers to run
// calculations in multiple threads, or possibly even use WebGL (via gpu.js?) to
// utilze graphics acceleration. in addition, due to lack of arbitrary-precision
// arithmetic, if you zoom in far enough, eventually rounding errors will cause
// features to become unresolvable.

var canvas = document.getElementById("canvas");
var ctx = canvas.getContext("2d");
ctx.fillStyle = "white";
ctx.textBaseline = "top";
ctx.textAlign = "left";
var [Nx, Ny] = [canvas.width, canvas.height];
var N = Nx * Ny;
var [cx, cy] = [0.5, 0.5];
var R2 = Math.pow(5., 2);
var TAU = 2*Math.PI;
var maxIter = 500;

// keeps track of the zoom level and central
// frame coordinates. also provides utility for
// converting from canvas pixel indices to frame
// coordinates
var zoomFrame = {
    zoom: 2./Nx,
    x: 0.,
    y: 0.,
    pixToCoord: function(i, j) {
        var x = this.zoom * (i - Nx/2.) + this.x;
        var y = this.zoom * (j - Nx/2.) + this.y;
        return [x, y];
    }
}

// keep track of where the mouse is pointing on the
// canvas so we can 1) update the coordinate display,
// and 2) provide a point to zoom into/out of when
// scrolling
var mouseCoords = {i: null, j: null};
var xCoord = document.getElementById("x-coord");
var yCoord = document.getElementById("y-coord");
canvas.onmousemove = function(e) {
    var mc = mouseCoords;
    var zf = zoomFrame;
    // flip the y axis coordinate so (0, 0) is bottom left
    [mc.i, mc.j] = [e.offsetX, Ny - e.offsetY];
    var [x, y] = zf.pixToCoord(mc.i, mc.j);
    xCoord.innerHTML = x.toExponential();
    yCoord.innerHTML = y.toExponential();
}

// updates frame and redraws image on scroll
canvas.onwheel = function(e) {
    var zf = zoomFrame;
    var mc = mouseCoords;
    var [oldX, oldY] = zf.pixToCoord(mc.i, mc.j);
    var newZoom = zf.zoom * Math.exp(0.1*e.deltaY);
    // this transformation ensures that the canvas pixel
    // the mouse is pointing at maps to the same (x,y)
    // coordinate in both the old/new frames for "smooth"
    // zooming
    zf.x = oldX - newZoom * (mc.i - Nx/2.);
    zf.y = oldY - newZoom * (mc.j - Ny/2.);
    zf.zoom = newZoom;
    draw()
}


// Sliders for chosing value of c
var cxSlider = document.getElementById("cx-slider");
cxSlider.value = cx;
cxSlider.onchange = function(e) {
    cx = parseFloat(cxSlider.value);
    draw();
}

var cySlider = document.getElementById("cy-slider");
cySlider.value = cy;
cySlider.onchange = function(e) {
    cy = parseFloat(cySlider.value);
    draw();
}


// calculates the "escape" iterations for any
// complex number z = x + i*y given the recurrence
// relation z <- z^2 + c
function mandelbrot(x, y, cx_, cy_) {
    var x_
    for (var i=0; i<maxIter; i++) {
        if (x*x + y*y > R2) {
            return i;
        } else {
            x_ = x;
            x = x*x - y*y + cx_;
            y = 2*x_*y + cy_;
        }
    }
    return 0;
}

// maps the escape iteration to an rgb color triplet
function cmap(i) {
    var t = TAU*i/20;
    return [
        // irrational frequencies ensures that no color is repeated
        Math.floor(255*Math.sin(t)), // red
        Math.floor(255*Math.sin(t*Math.SQRT2)), // green
        Math.floor(255*Math.sin(t*Math.PI)) // blue
    ];
}

// redraws the fractal on the canvas
function draw() {
    var start = Date.now();
    var img = new ImageData(Nx, Ny);
    var offset, x, y;
    for (var i=0; i<Nx; i++) {
        for (var j=0; j<Ny; j++) {
            offset = 4*(Ny*(Ny-j-1)+i);
            [x, y] = zoomFrame.pixToCoord(i, j);
            var z = mandelbrot(x, y, cx, cy)
            var [r, g, b] = cmap(z);
            img.data[offset] = Math.floor(r);
            img.data[offset+1] = Math.floor(g);
            img.data[offset+2] = Math.floor(b);
            img.data[offset+3] = 255;
        }
    }
    ctx.putImageData(img, 0, 0);
    var z = Math.log10(2./(Nx*zoomFrame.zoom)).toFixed(2);
    ctx.fillText(`Zoom: 10^${z}x`, 0.02*canvas.width, 0.02*canvas.height);
    var renderTime = Date.now() - start;
    ctx.fillText(`Render time: ${renderTime} ms`, 0.02*canvas.width, 0.02*canvas.height + 15);
}

// draw first frame
draw();