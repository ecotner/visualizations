// initialize canvases
var L = 100;
var B = 1e-3;
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
    var dt = 1e-5;
    var N = 100;
    var dx = 0.1;
    var wf = new Module.WaveFunction(
        L, L, dx, dx*L/2., dx*L/4., 5.0/(dx*L), 5.0/(dx*L), 1., "dirichlet");
    wf.stepFTCS(dt, B);
    var t1 = 0.
    function asdf(t2) {
        draw(wf);
        wf.multiStep(N, dt, B);
        t_frame = t2 - t1;
        t1 = t2;
        console.log(t_frame);
        // console.log(wf.abs2(dx*L/2, dx*L/2));
        window.requestAnimationFrame(asdf);
    }
    window.requestAnimationFrame(function (t2) {
        t1 = t2;
        window.requestAnimationFrame(asdf);
    })
}

var Module = {onRuntimeInitialized: run};