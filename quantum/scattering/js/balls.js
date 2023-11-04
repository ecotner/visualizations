import { CanvasManager } from "./canvas.js";

function doesBallOverlapWall(ball, wall) {
    [xb, yb] = [ball[0], ball[1]]
    ;
}

export class BallGame {
    /**
     * 
     * @param {number} Nx The number of grid points in the x direction
     * @param {number} Ny The number of grid points in the y direction
     * @param {number} r The radius of the balls in grid units
     */
    constructor(Nx, Ny, r) {
        this.Nx = Nx
        this.Ny = Ny
        this.r = r
        this.cm = new CanvasManager(Ny, Nx)
        this.ballStates = []
        this.wallShapes = []
        this.prev_time = null
        this.mousedownCoords = null

        // add event listener to create balls when canvas is clicked
        // give them initial velocity by holding mouse down and moving
        this.cm.canvas.addEventListener("mousedown", (event) => {
            this.mousedownCoords = [event.offsetX, this.cm.height - event.offsetY]
        })
        this.cm.canvas.addEventListener("mouseup", (event) => {
            const [x1, y1] = this.mousedownCoords
            const [x2, y2] = [event.offsetX, this.cm.height - event.offsetY]
            const [vx, vy] = [x2 - x1, y2 - y1] // velocity units are px/sec
            this.addBall(x1, y1, vx, vy)
        })
    }

    addBall(x, y, vx, vy) {
        this.ballStates.push([x, y, vx, vy])
        console.log(`added ball at (${x}, ${y}); there are ${this.ballStates.length} balls`)
    }

    addWall(x1, y1, x2, y2) {
        // swap x and y coordinates if not in correct format
        if (x1 > x2) {
            [x1, x2] = [x2, x1]
        }
        if (y1 > y2) {
            [y1, y2] = [y2, y1]
        }
        // make sure wall isn't being created inside a pre-existing ball
        this.wallShapes.push([x1, y1, x2, y2])
        console.log(`added wall with corners (${x1}, ${y1}) and (${x2}, ${y2})`)
    }

    render() {
        const ctx = this.cm.ctx
        ctx.clearRect(0, 0, this.Nx, this.Ny)
        for (let ball of this.ballStates) {
            ctx.beginPath()
            ctx.arc(ball[0], ball[1], this.r, 0, 2 * Math.PI)
            ctx.fill()
        }
    }
    
    updateBallStates(dt) {
        const [ax, ay, r] = [0, -10, this.r] // acceleration is in px/sec^2
        // free ballistic motion update
        for (let i = 0; i < this.ballStates.length; i++) {
            const ball = this.ballStates[i]
            ball[0] = ball[0] + ball[2]*dt + 0.5*ax*dt*dt // x
            ball[1] = ball[1] + ball[3]*dt + 0.5*ay*dt*dt // y
            ball[2] = ball[2] + ax*dt // vx
            ball[3] = ball[3] + ay*dt // vy
        }
        // delete balls that have fallen off screen
        for (let i = this.ballStates.length - 1; i >= 0; i--) {
            const ball = this.ballStates[i]
            if (ball[1] + r < 0 || ball[0] + r < 0 || ball[0] - r > this.Nx || ball[1] - r > this.Ny) {
                this.ballStates.splice(i, 1)
            }
        }
        // look for and correct collision trajectories between balls
        for (let i = 0; i< this.ballStates.length; i++) {
            // undo trajectory to the point where the balls are touching exactly; figure out how long that takes (dt = dt1 + dt2)
            // find new velocities using conservation of energy/momentum
            // evolve trajectory post-collision (for dt2)
        }
    }

    run() {
        const frameIdx = window.requestAnimationFrame((time) => {
            // console.log(`rendering frame ${frameIdx} at time ${time}`)
            if (this.prev_time === null) this.prev_time = time
            let dt = (time - this.prev_time) / 1000.0 // time elapsed in seconds
            this.prev_time = time
            this.updateBallStates(dt)
            this.render()
            this.run()
        })
    }
}