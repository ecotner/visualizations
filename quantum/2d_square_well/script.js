// initialize canvases
var L = 50;
var qB = 0;
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

function draw(wf) {
    memCtx.clearRect(0, 0, canvasWidth, canvasHeight);
    var id = memCtx.getImageData(0, 0, canvasWidth, canvasHeight);
    var pixels = id.data;
    for (var i=0; i<canvasWidth; i++) {
        for (var j=0; j<canvasHeight; j++) {
            var offset = 4*(id.width * j + i);
            var psi = wf.abs2(i * scaleX, j * scaleY);
            // var intensity = 1e1*(Math.log10(psi) + 15);
            var intensity = Math.floor(1000*psi);
            // pixels[offset] = intensity;
            pixels[offset] = 255 - intensity;
            pixels[offset+1] = 255 - intensity;
            pixels[offset+2] = 255 - intensity;
            pixels[offset+3] = 255;
        }
    }
    memCtx.putImageData(id, 0, 0);
    ctx.drawImage(memCanvas, 0, 0);
}

function run() {
    var dt = 5e-1;
    var N = 100;
    var m = 10;
    var dx = 1.;
    var wf = new Module.WaveFunction(
        L, L, dx, L/2., L/4., 5./L, 5./L, L/10., "dirichlet");
    wf.stepFTCS(dt, m, qB);
    function asdf() {
        draw(wf);
        wf.multiStep(N, dt, m, qB);
        // console.log(wf.abs2(L/2, L/2));
        window.requestAnimationFrame(asdf);
    }
    asdf();
}

var Module = {onRuntimeInitialized: run};