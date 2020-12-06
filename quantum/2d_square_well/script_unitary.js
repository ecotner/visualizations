var N = 100;
var dx = 0.1;
var L = (N-1)*dx;
var dt = 1e-3;
var B = 1e-1;
var x0 = 0.5*L;
var y0 = 0.5*L;
var s0 = 1.;
var tau = 2*Math.PI;
var px0 = 1*tau/L;
var py0 = 1*tau/L;


// initialize canvases
var canvas = document.getElementById("canvas");
var ctx = canvas.getContext("2d");
var canvasWidth = canvas.width;
var canvasHeight = canvas.height;
var scaleX = L / canvasWidth;
var scaleY = L / canvasHeight;
var memCanvas = document.createElement("canvas");
Object.assign(memCanvas, {width: canvasWidth, height: canvasHeight})
var memCtx = memCanvas.getContext("2d");
ctx.translate(0, canvasHeight);
ctx.scale(1, -1)

var MAX_PSI2 = 1;
function draw(wf) {
    memCtx.clearRect(0, 0, canvasWidth, canvasHeight);
    var id = memCtx.getImageData(0, 0, canvasWidth, canvasHeight);
    var pixels = id.data;
    var _max_psi2 = 0;
    for (var i=0; i<canvasWidth; i++) {
        for (var j=0; j<canvasHeight; j++) {
            var offset = 4*(id.width * j + i);
            var psi2 = wf.abs2(i * scaleX, j * scaleY);
            _max_psi2 = Math.max(psi2, _max_psi2);
            var intensity = Math.floor(255*psi2/MAX_PSI2);
            pixels[offset+1] = intensity;
            // pixels[offset] = 255 - intensity;
            // pixels[offset+1] = 255 - intensity;
            // pixels[offset+2] = 255 - intensity;
            pixels[offset+3] = 255;
        }
    }
    MAX_PSI2 = _max_psi2;
    console.log(MAX_PSI2);
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
}

var wf;
function init() {
    wf = new Module.WaveFunction(N, N, dx, dt, B, x0, y0, s0, px0, py0);
    FRAME_ID = startAnimation(wf);
}

var Module = {onRuntimeInitialized: init};

// buttons
var startBtn = document.getElementById("start-button");
startBtn.onclick = function(event) {
    startAnimation(wf)
}

var pauseBtn = document.getElementById("pause-button");
pauseBtn.onclick = function(event) {
    window.cancelAnimationFrame(FRAME_ID);
}

var resetBtn = document.getElementById("reset-button");
resetBtn.onclick = function(event) {
    window.cancelAnimationFrame(FRAME_ID);
    wf.delete();
    wf = new Module.WaveFunction(N, N, dx, dt, B, x0, y0, s0, px0, py0);
    draw(wf);
}