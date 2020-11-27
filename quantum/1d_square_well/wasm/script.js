var canvas = document.getElementById("canvas");
var ctx = canvas.getContext("2d");
ctx.strokeStyle = "red";
ctx.translate(0, canvas.height);
ctx.scale(1, -1);

function draw(wf) {
    var L = wf.getL();
    var x
    var yscale = 0.6 * canvas.height;
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.beginPath();
    ctx.moveTo(0, 0.02*yscale * (Math.log(wf.getPsiAbsSq(0)) + 50));
    
    for (var i=1; i<canvas.width; i++) {
        x = (i * L)/(canvas.width-1);
        ctx.lineTo(i, 0.02*yscale * (Math.log(wf.getPsiAbsSq(x)) + 50));
    }
    ctx.stroke();
    wf.multiStep(500);
    window.requestAnimationFrame(function () {draw(wf)});
}

function run() {
    var wf = new Module.WaveFunction(500, 1e-3, 5.0);
    draw(wf);
}

var Module = {onRuntimeInitialized: run};