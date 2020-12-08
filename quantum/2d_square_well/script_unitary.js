// declare parameters in global scope so they can be modified
var N = 100;
var dx = 0.1;
var L = (N-1)*dx;
var dt = 1e-3;
var B = 0.1;
var E = 0;
var thetaE = 0;
var x0 = 0.5;
var y0 = 0.5;
var s0 = 0.1;
var tau = 2*Math.PI;
var k = tau/L;
var p0 = 0;
var theta0 = 45;
var simState = "paused";

// initialize canvases
var canvas = document.getElementById("canvas"); // canvas to render
var ctx = canvas.getContext("2d");
var memCanvas = document.createElement("canvas"); // in-memory canvas only
Object.assign(memCanvas, {width: canvas.width, height: canvas.height})
var memCtx = memCanvas.getContext("2d");
ctx.translate(0, canvas.height);
ctx.scale(1, -1)

var MAX_PSI2 = 1; // used for rescaling pixel intensities
// function for drawing wave function on canvas
function draw(wf) {
    memCtx.clearRect(0, 0, canvas.width, canvas.height);
    var id = memCtx.getImageData(0, 0, canvas.width, canvas.height);
    var pixels = id.data;
    var _max_psi2 = 0;
    var scaleX = L / canvas.width;
    var scaleY = L / canvas.height;
    for (var i=0; i<canvas.width; i++) {
        for (var j=0; j<canvas.height; j++) {
            var offset = 4*(id.width * j + i);
            var psi2 = wf.abs2(i * scaleX, j * scaleY);
            _max_psi2 = Math.max(psi2, _max_psi2);
            var intensity = Math.floor(256*psi2/MAX_PSI2);
            pixels[offset+1] = intensity;
            // pixels[offset] = 255 - intensity;
            // pixels[offset+1] = 255 - intensity;
            // pixels[offset+2] = 255 - intensity;
            pixels[offset+3] = 255;
        }
    }
    MAX_PSI2 = _max_psi2;
    memCtx.putImageData(id, 0, 0);
    ctx.drawImage(memCanvas, 0, 0);
}

var FRAME_ID;
function startAnimation(wf) {
    draw(wf);
    wf.step(5);
    FRAME_ID = window.requestAnimationFrame(function() {
        startAnimation(wf);
    });
    console.log(wf.probMass());
}

// parameter controls
var BSlider = document.getElementById("B-slider");
BSlider.onchange = function(event) {
    B = parseFloat(BSlider.value);
    var text = B.toFixed(2);
    if (B > 0) text = ("+" + text).padStart(5, " ")
    document.getElementById("B-value").innerHTML = text;
}

var NSlider = document.getElementById("N-slider");
NSlider.onchange = function(event) {
    N = parseInt(NSlider.value);
    L = dx*(N-1);
    var text = N.toString().padStart(3, " ");
    text += " (L: " + L.toFixed(2).padStart(5, " ") + "/m)";
    document.getElementById("N-value").innerHTML = text;
}

var dtSlider = document.getElementById("dt-slider");
dtSlider.onchange = function(event) {
    dt = Math.pow(10, parseFloat(dtSlider.value));
    var text = dt.toExponential(2).padStart(6, " ");
    document.getElementById("dt-value").innerHTML = text;
}

var dxSlider = document.getElementById("dx-slider");
dxSlider.onchange = function(event) {
    dx = Math.pow(10, parseFloat(dxSlider.value));
    L = dx*(N-1);
    var text = dx.toExponential(2).padStart(6, " ");
    document.getElementById("dx-value").innerHTML = text;
}

var x0Slider = document.getElementById("x0-slider");
x0Slider.onchange = function(event) {
    x0 = parseFloat(x0Slider.value);
    var text = x0.toFixed(2).padStart(4, " ");
    document.getElementById("x0-value").innerHTML = text;
}

var y0Slider = document.getElementById("y0-slider");
y0Slider.onchange = function(event) {
    y0 = parseFloat(y0Slider.value);
    var text = y0.toFixed(2).padStart(4, " ");
    document.getElementById("y0-value").innerHTML = text;
}

var s0Slider = document.getElementById("s0-slider");
s0Slider.onchange = function(event) {
    s0 = parseFloat(s0Slider.value);
    var text = s0.toFixed(2).padStart(4, " ");
    document.getElementById("s0-value").innerHTML = text;
}

var p0Slider = document.getElementById("p0-slider");
p0Slider.onchange = function(event) {
    p0 = parseFloat(p0Slider.value);
    var text = p0.toFixed(2).padStart(5, " ");
    document.getElementById("p0-value").innerHTML = text;
}

var theta0Slider = document.getElementById("theta0-slider");
theta0Slider.onchange = function(event) {
    theta0 = parseInt(theta0Slider.value);
    var text = theta0.toString().padStart(3, " ");
    document.getElementById("theta0-value").innerHTML = text;
}

var ESlider = document.getElementById("E-slider");
ESlider.onchange = function(event) {
    E = parseFloat(ESlider.value);
    var text = E.toFixed(2).padStart(5, " ");
    document.getElementById("E-value").innerHTML = text;
}

var EthetaSlider = document.getElementById("E-theta-slider");
EthetaSlider.onchange = function(event) {
    thetaE = parseInt(EthetaSlider.value);
    var text = thetaE.toString().padStart(3, " ");
    document.getElementById("E-theta-value").innerHTML = text;
}

var wf; // wave function object exists in global scope so other functions can modify it
function init() {
    // prepare sliders
    for (var vars of [
        [BSlider, B], [NSlider, N], [dtSlider, Math.log10(dt)], [dxSlider, Math.log10(dx)],
        [x0Slider, x0], [y0Slider, y0], [s0Slider, s0], [p0Slider, p0], [theta0Slider, theta0],
        [ESlider, E], [EthetaSlider, thetaE]
    ]) {
        [slider, x] = vars;
        slider.value = x;
        slider.onchange();
    }
    // initialize simulation and start animation
    var px = p0 * k * Math.cos(theta0 * Math.PI / 180);
    var py = p0 * k * Math.sin(theta0 * Math.PI / 180);
    var Ex = E * Math.cos(thetaE);
    var Ey = E * Math.sin(thetaE);
    wf = new Module.WaveFunction(N, N, dx, dt, B, Ex, Ey, x0*L, y0*L, s0*L, px, -py);
    FRAME_ID = startAnimation(wf);
}

var Module = {onRuntimeInitialized: init};

// control panel buttons
var startBtn = document.getElementById("start-button");
startBtn.onclick = function(event) {
    simState = "running";
    startAnimation(wf)
}

var pauseBtn = document.getElementById("pause-button");
pauseBtn.onclick = function(event) {
    simState = "paused";
    window.cancelAnimationFrame(FRAME_ID);
}

var resetBtn = document.getElementById("reset-button");
resetBtn.onclick = function(event) {
    simState = "paused";
    window.cancelAnimationFrame(FRAME_ID);
    wf.delete();
    var px = p0 * k * Math.cos(theta0 * Math.PI / 180);
    var py = p0 * k * Math.sin(theta0 * Math.PI / 180);
    var Ex = E * Math.cos(thetaE);
    var Ey = E * Math.sin(thetaE);
    wf = new Module.WaveFunction(N, N, dx, dt, B, Ex, Ey, x0*L, y0*L, s0*L, px, -py);
    draw(wf);
}