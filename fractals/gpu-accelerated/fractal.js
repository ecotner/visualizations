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
const gpu = new GPU();
var [Nx, Ny] = [canvas.width, canvas.height];
var N = Nx * Ny;
var [cx, cy] = [0.5, 0.5];
var R2 = Math.pow(5., 2);
var TAU = 2*Math.PI;
var maxIter = 500.;

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
    xCoord.innerHTML = x.toExponential(3);
    yCoord.innerHTML = y.toExponential(3);
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


// calculates the "escape" iterations for any
// complex number z = x + i*y given the recurrence
// relation z <- z^2 + c
function mandelbrot(x, y) {
    for (var i=0; i<maxIter; i++) {
        if (x*x + y*y > R2) {
            return i;
        } else {
            [x, y] = [x*x - y*y + cx, 2*x*y + cy];
        }
    }
    return 0;
}

const mandelbrotGPU = gpu.createKernel(function(x0, y0, dx, dy) {
    var x = x0 + dx * this.thread.x;
    var y = y0 + dy * this.thread.y;
    var xTemp = 0;
    var ans = 0, cx = this.constants.cx, cy = this.constants.cy;
    var R2 = this.constants.R2;
    for (var i=0; i<this.constants.maxIter; i++) {
        if (x*x + y*y > R2) {
            ans = i;
            break;
        } else {
            xTemp = x;
            x = x*x - y*y + cx;
            y = 2*xTemp*y + cy;
        }
    }
    return ans;
}, {
    constants: {
        maxIter: maxIter,
        R2: R2,
        cx: cx,
        cy: cy
    },
    loopMaxIterations: maxIter+1,
}).setOutput([Nx, Ny]);

// maps the escape iteration to an rgb color triplet
function cmap(i) {
    var t = TAU*i/20;
    return [
        // irrational frequencies ensures that no color is repeated
        255*Math.sin(t), // red
        255*Math.sin(t*Math.SQRT2), // green
        255*Math.sin(t*Math.PI) // blue
    ];
}

// redraws the fractal on the canvas
function draw() {
    var start = Date.now();
    var img = new ImageData(Nx, Ny);
    var offset, x, y, x0, y0;
    [x0, y0] = zoomFrame.pixToCoord(0, 0);
    var dx = zoomFrame.zoom
    var vals = mandelbrotGPU(x0, y0, dx, dx);
    for (var i=0; i<Nx; i++) {
        for (var j=0; j<Ny; j++) {
            offset = 4*(Ny*(Ny-j-1)+i);
            var [r, g, b] = cmap(vals[j][i]);
            img.data[offset] = Math.floor(r);
            img.data[offset+1] = Math.floor(g);
            img.data[offset+2] = Math.floor(b);
            img.data[offset+3] = 255;
        }
    }
    ctx.putImageData(img, 0, 0);
    var z = (2./(Nx*zoomFrame.zoom)).toFixed(1);
    ctx.fillText(`Zoom: ${z}x`, 0.02*canvas.width, 0.02*canvas.height);
    var renderTime = Date.now() - start;
    ctx.fillText(`Render time: ${renderTime} ms`, 0.02*canvas.width, 0.02*canvas.height + 15);
}

// draw first frame
draw();