import { CanvasManager } from "./canvas.js";

export class QuantumGame {
    constructor(Nx, Ny, m) {
        this.Nx = Nx
        this.Ny = Ny
        this.m = m
        this.cm = new CanvasManager(Ny, Nx)
        this.prev_time = null
        this.mousedownCoords = null
        this.objCreationState = "balls" // can be either "balls" or "walls"
        this.maxProbDensityExpWeight = 0 // exponentially-weighted max probability density
        this.alphaExpWeight = 0.999 // exponential weighting factor

        // wavefunction components
        this.A = [] // real part
        this.B = [] // imaginary part
        // states at previous time; required for using centered-time evolution scheme because
        // forward-time euler is not stable
        this.A_prev = []
        this.B_prev = []
        // memory buffers for writing next state without overwriting current/previous states
        this.A_next = []
        this.B_next = []
        // flattened rasterization of |psi|^2
        this.rasterNorm2 = []
        // fill with all zeros
        for (let i=0; i < Nx; i++) {
            this.A.push([])
            this.B.push([])
            this.A_prev.push([])
            this.B_prev.push([])
            this.A_next.push([])
            this.B_next.push([])
            for (let j=0; j < Ny; j++) {
                this.A[i].push(0.0)
                this.B[i].push(0.0)
                this.A_prev[i].push(0.0)
                this.B_prev[i].push(0.0)
                this.A_next[i].push(0.0)
                this.B_next[i].push(0.0)
                this.rasterNorm2.push(0.0)
            }
        }

        // add event listener to create wavepackets/walls when canvas is clicked
        // give wavepackets initial velocity by holding mouse down and moving
        this.cm.canvas.addEventListener("mousedown", (event) => {
            this.mousedownCoords = [event.offsetX, this.cm.height - event.offsetY]
        })
        this.cm.canvas.addEventListener("mouseup", (event) => {
            const [x1, y1] = this.mousedownCoords
            const [x2, y2] = [event.offsetX, this.cm.height - event.offsetY]
            const [vx, vy] = [x2 - x1, y2 - y1] // velocity units are px/sec
            if (this.objCreationState === "balls") this.addWavePacket(x1, y1, vx, vy)
            else if (this.objCreationState === "walls") this.addWall(x1, y1, x2, y2)
            // else console.log(`unrecognized object creation state '${this.objCreationState}'`)
        })

        // this.addWavePacket(Nx/2, Ny/2, 2, 0)
    }

    addWavePacket(xp, yp, vx, vy) {
        const [Nx, Ny, A, B, A_prev, B_prev, m] = [
            this.Nx, this.Ny, this.A, this.B, this.A_prev, this.B_prev, this.m
        ]
        const [px, py] = [m * vx, m * vy]
        const width = 1/m
        const expCoeff = 1/(2*width*width)
        const normCoeff = 1/(width*Math.sqrt(2*Math.PI))
        // add prob. amplitude to previous time step
        for (let x=1; x < Nx-1; x++) {
            for (let y=1; y < Ny-1; y++) {
                const norm = normCoeff * Math.exp(-expCoeff*((x-xp)*(x-xp) + (y-yp)*(y-yp)))
                const phase = (px * x) + (py * y)
                A_prev[x][y] = A_prev[x][y] + norm * Math.cos(phase)
                B_prev[x][y] = B_prev[x][y] + norm * Math.sin(phase)
            }
        }
        // perform single forward euler step to make current state consistent
        const dt_2mdx2 = 0.010 / (2*this.m)
        for (let x=1; x < Nx-1; x++) {
            for (let y=1; y < Ny-1; y++) {
                A[x][y] = A_prev[x][y] - dt_2mdx2 * (B_prev[x-1][y] + B_prev[x+1][y] + B_prev[x][y-1] + B_prev[x][y+1] - 4*B_prev[x][y] )
                B[x][y] = B_prev[x][y] + dt_2mdx2 * (A_prev[x-1][y] + A_prev[x+1][y] + A_prev[x][y-1] + A_prev[x][y+1] - 4*A_prev[x][y] )
            }
        }
        console.log(`added wave packet at (${xp}, ${yp}) with velocity (${vx}, ${vy})`)
    }

    updateWavefunction(dt) {
        const [A, B, A_prev, B_prev, A_next, B_next, m, Nx, Ny] = [
            this.A, this.B, this.A_prev, this.B_prev, this.A_next, this.B_next,
            this.m, this.Nx, this.Ny
        ]
        const dt_2mdx2 = dt / (2*m)
        let probMass = 0.0
        // update wavefunction at grid points and measure total probability mass
        // iterate only over grid points [1, N-2] instead of [0, N-1] because boundaries are fixed to zero
        for (let x=1; x < Nx-1; x++) {
            for (let y=1; y < Ny-1; y++) {
                const An = A_prev[x][y] - dt_2mdx2 * (B[x-1][y] + B[x+1][y] + B[x][y-1] + B[x][y+1] - 4*B[x][y])
                const Bn = B_prev[x][y] + dt_2mdx2 * (A[x-1][y] + A[x+1][y] + A[x][y-1] + A[x][y+1] - 4*A[x][y])
                probMass += An*An + Bn*Bn
                A_next[x][y] = An
                B_next[x][y] = Bn
            }
        }
        // normalize the total probability mass to 1
        const norm = Math.sqrt(probMass)
        if (norm > 0) {
            for (let x=0; x < Nx; x++) {
                for (let y=0; y < Ny; y++) {
                    A_next[x][y] /= norm
                    B_next[x][y] /= norm
                }
            }
        }
        // swap around the array assignments
        this.A = A_next
        this.A_prev = A
        this.A_next = A_prev
        this.B = B_next
        this.B_prev = B
        this.B_next = B_prev
    }

    render() {
        const [A, B, Nx, Ny] = [this.A, this.B, this.Nx, this.Ny]
        // wipe the slate clean
        this.cm.ctx.clearRect(0, 0, this.Nx, this.Ny)
        // flatten the wavefunction arrays and compute norm
        const arr = this.rasterNorm2
        let probMass = 0.0
        let maxProbDensity = -1
        for (let x=0; x < Nx; x++) {
            for (let y=0; y < Ny; y++) {
                const i = Nx*(Ny-1-y) + x
                arr[i] = A[x][y]*A[x][y] + B[x][y]*B[x][y]
                probMass += arr[i]
                maxProbDensity = Math.max(arr[i], maxProbDensity)
            }
        }
        // adjust the rasterized image pixel brightness scale so that the point of highest
        // probability density (exponentiall weighted over time to smooth things out) is
        // scaled to 1
        if (probMass > 0) {
            const [alpha, weight] = [this.alphaExpWeight, this.maxProbDensityExpWeight]
            this.maxProbDensityExpWeight = alpha * weight + (1 - alpha) * maxProbDensity
            for (let i=0; i < arr.length; i++) {
                arr[i] = arr[i] / this.maxProbDensityExpWeight
            }
        }
        this.cm.render_raster_float(arr)
        // console.log(`total prob. mass: ${probMass}`)
    }

    run() {
        const frameIdx = window.requestAnimationFrame((time) => {
            // console.log(`rendering frame ${frameIdx} at time ${time}`)
            if (this.prev_time === null) this.prev_time = time
            const dt = (time - this.prev_time) / 1000.0 // time elapsed in seconds
            this.prev_time = time
            const N = 10
            for (let i=0; i<N; i++) {
                this.updateWavefunction(dt/N)
            }
            this.render()
            this.run()
        })
    }
}