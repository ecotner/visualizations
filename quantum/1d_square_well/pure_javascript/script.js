// define mathematical constants for initializing wave function
const L = 300;
const tau = 2*Math.PI;
const sigma = L/10;
const p = 5*tau/L;

// initialize gpu.js instance
const gpu = new GPU();

// get canvas, switch to euclidean coordinates
const canvas = document.getElementById("canvas");
const ctx = canvas.getContext("2d");
ctx.translate(0, canvas.height);
ctx.scale(1, -1);

// define some logical comparison functions based on floating-point math
function ge(x, y) {
    return Math.sign(Math.sign(x - y) + 1);
}
function le(x, y) {
    return Math.sign(Math.sign(y - x) + 1);
}
gpu.addFunction(ge);
gpu.addFunction(le);

// defines the time evolution single-step update for the wave function
const update = gpu.createKernel(function(psi, psi_prev) {
    var [x, y] = [this.thread.x, this.thread.y];
    var psi_next = psi_prev[y][x] - (
        Math.sign(y-0.5) * (this.constants.dt/2)
            * (psi[1-y][x+1] - 2*psi[1-y][x] + psi[1-y][x-1])
    )
    // conditional branches degrade performance on GPU, so I will
    // apply boundary conditions by masking
    var mask = le(x, this.constants.L-1) * ge(x, 1);
    return mask*psi_next;
}, {
    constants: {L: L, dt: 1.0},
    pipeline: true, // allows kernels to pass data without sending to cpu
    immutable: true, // manual memory management; needs to be enabled?
    output: [L+1, 2], // shape of output
})

// fuses individual updates into multiple
// much faster, but also less numerically stable for some reason?
const multiUpdate = gpu.combineKernels(update, function(psi, psi_prev) {
    for (var i=0; i<100; i++) {
        var psi_next = update(psi, psi_prev);
        if (i>1) psi_prev.delete();
        [psi, psi_prev] = [psi_next, psi];
    }
    return [psi_next, psi];
})

// initialize the wave function real/imaginary parts
function initPsi() {
    let psi = [[],[]]
    const N = Math.sqrt(0.7*300);
    for (var x=0; x<L+1; x++) {
        rho = N * Math.exp(-Math.pow(x-L/2, 2)/(2*Math.pow(sigma, 2)));
        phase = [Math.cos(-p*x), Math.sin(-p*x)];
        for (var i=0; i<2; i++) {
            psi[i].push(rho*phase[i]);
        }
    }
    for (var i=0; i<2; i++){
        psi[i][0] = 0; psi[i][L+1] = 0;
    }
    return psi;
}

// draw the absolute square of the wave function (probability
// density) on the canvas
function draw(psi) {
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    var psi_re = psi[0];
    var psi_im = psi[1];
    ctx.beginPath();
    ctx.moveTo(0, 0);
    // iterate over pixels along width of canvas
    var y_re = 0;
    var y_im = 0;
    for (var i=1; i<=canvas.width; i++) {
        // find closest bounding grid points
        var x = L * i / canvas.width;
        var [x1, x2] = [Math.floor(x), Math.ceil(x)];
        if (x1 == x2) { // falls exactly on a grid point
            y_re = psi_re[x1]
            y_im = psi_im[x1];
        } else { // interpolate between grid points
            var [y1, y2] = [psi_re[x1], psi_re[x2]];
            y_re = y1 + (y2-y1)/(x2-x1)*(x-x1);
            [y1, y2] = [psi_im[x1], psi_im[x2]];
            y_im = y1 + (y2-y1)/(x2-x1)*(x-x1);
        }
        // draw line to next pixel in canvas
        var rho2 = y_re*y_re + y_im*y_im;
        ctx.lineTo(i, rho2);
    }
    ctx.stroke();
}

var psi = initPsi();
var psi_prev = initPsi();
var psi_next;
var i = 0;
function run() {
    let tic = new Date();
    for (var j=0; j<100; j++) {
        var psi_next = update(psi, psi_prev);
        // delete old textures as they are no longer needed to avoid memory leaks
        if (i>1) psi_prev.delete();
        [psi, psi_prev] = [psi_next, psi];
        i++;
    }
    let toc = new Date();
    console.log(toc - tic) // observe the frame rate in the console
    draw(psi.toArray()); // draw probability density
    window.requestAnimationFrame(run)
}
run()
// console.log(psi_next.toArray());
// [psi_next, psi, psi_prev].map(e => e.delete());