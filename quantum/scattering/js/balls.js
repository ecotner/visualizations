import { CanvasManager } from "./canvas.js";

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
        this.ax = 0
        this.ay = -100
        this.objCreationState = "balls" // can be either "balls" or "walls"

        // add event listener to create balls when canvas is clicked
        // give them initial velocity by holding mouse down and moving
        this.cm.canvas.addEventListener("mousedown", (event) => {
            this.mousedownCoords = [event.offsetX, this.cm.height - event.offsetY]
        })
        this.cm.canvas.addEventListener("mouseup", (event) => {
            const [x1, y1] = this.mousedownCoords
            const [x2, y2] = [event.offsetX, this.cm.height - event.offsetY]
            const [vx, vy] = [x2 - x1, y2 - y1] // velocity units are px/sec
            if (this.objCreationState === "balls") this.addBall(x1, y1, vx, vy)
            else if (this.objCreationState === "walls") this.addWall(x1, y1, x2, y2)
            else console.log(`unrecognized object creation state '${this.objCreationState}'`)
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
        // wipe the slate clean
        ctx.clearRect(0, 0, this.Nx, this.Ny)
        // draw balls
        for (let ball of this.ballStates) {
            ctx.beginPath()
            ctx.arc(ball[0], ball[1], this.r, 0, 2 * Math.PI)
            ctx.fill()
        }
        // draw walls
        for (let wall of this.wallShapes) {
            ctx.beginPath()
            ctx.strokeRect(wall[0], wall[1], wall[2] - wall[0], wall[3] - wall[1])
        }
    }

    updateBallisticTrajectory(dt) {
        const [ax, ay] = [this.ax, this.ay]
        // free ballistic motion update
        for (let i = 0; i < this.ballStates.length; i++) {
            const ball = this.ballStates[i]
            ball[0] = ball[0] + ball[2]*dt + 0.5*ax*dt*dt // x
            ball[1] = ball[1] + ball[3]*dt + 0.5*ay*dt*dt // y
            ball[2] = ball[2] + ax*dt // vx
            ball[3] = ball[3] + ay*dt // vy
        }
    }

    deleteOffScreenBalls() {
        // delete balls that have fallen off screen
        const r = this.r
        for (let i = this.ballStates.length - 1; i >= 0; i--) {
            const ball = this.ballStates[i]
            if (ball[1] + r < 0 || ball[0] + r < 0 || ball[0] - r > this.Nx || ball[1] - r > this.Ny) {
                this.ballStates.splice(i, 1)
            }
        }
    }

    handleBallBallCollisions(dt) {
        const [ax, ay, r] = [this.ax, this.ay, this.r]
        // look for and correct collision trajectories between balls
        for (let i = 0; i < this.ballStates.length; i++) {
            for (let j = i+1; j < this.ballStates.length; j++) {
                const [b1, b2] = [this.ballStates[i], this.ballStates[j]]
                let [dx, dy] = [b1[0] - b2[0], b1[1] - b2[1]]
                if (dx*dx + dy*dy > 4*r*r) continue // skip if balls are not touching
                // undo trajectories to the previous time step
                // ball 1
                b1[3] = b1[3] - ay*dt // vy
                b1[2] = b1[2] - ax*dt // vx
                b1[1] = b1[1] - b1[3]*dt - 0.5*ay*dt*dt // y
                b1[0] = b1[0] - b1[2]*dt - 0.5*ax*dt*dt // x
                // ball 2
                b2[3] = b2[3] - ay*dt // vy
                b2[2] = b2[2] - ax*dt // vx
                b2[1] = b2[1] - b2[3]*dt - 0.5*ay*dt*dt // y
                b2[0] = b2[0] - b2[2]*dt - 0.5*ax*dt*dt // x
                ; // need this here for some reason
                // run trajectory forward to the point where the balls are touching exactly; figure out how long that takes (dt = dt1 + dt2)
                [dx, dy] = [b1[0] - b2[0], b1[1] - b2[1]]
                const [dvx, dvy] = [b1[2] - b2[2], b1[3] - b2[3]]
                const [a, b, c] = [dvx*dvx + dvy*dvy, 2*(dx*dvx + dy*dvy), dx*dx + dy*dy - 4*r*r]
                const dt1 = (-b - Math.sqrt(b*b - 4*a*c)) / (2*a)
                // ball 1
                b1[0] = b1[0] + b1[2]*dt1 + 0.5*ax*dt1*dt1 // x
                b1[1] = b1[1] + b1[3]*dt1 + 0.5*ay*dt1*dt1 // y
                b1[2] = b1[2] + ax*dt1 // vx
                b1[3] = b1[3] + ay*dt1 // vy
                // ball 2
                b2[0] = b2[0] + b2[2]*dt1 + 0.5*ax*dt1*dt1 // x
                b2[1] = b2[1] + b2[3]*dt1 + 0.5*ay*dt1*dt1 // y
                b2[2] = b2[2] + ax*dt1 // vx
                b2[3] = b2[3] + ay*dt1 // vy
                // find new post-collision velocities using conservation of energy and momentum (linear + angular)
                // see https://en.wikipedia.org/wiki/Elastic_collision#Two-dimensional_collision_with_two_moving_objects
                // and https://williamecraver.wixsite.com/elastic-equations
                const [x1, y1, x2, y2, vx1, vy1, vx2, vy2] = [b1[0], b1[1], b2[0], b2[1], b1[2], b1[3], b2[2], b2[3]]
                const [v1, v2, theta1, theta2, phi] = [
                    Math.sqrt(vx1*vx1 + vy1*vy1),
                    Math.sqrt(vx2*vx2 + vy2*vy2),
                    Math.atan2(vy1, vx1),
                    Math.atan2(vy2, vx2),
                    Math.atan2((y2 - y1), x2 - x1)
                ]
                b1[2] = v2*Math.cos(theta2-phi)*Math.cos(phi) + v1*Math.sin(theta1-phi)*Math.cos(phi+0.5*Math.PI)
                b1[3] = v2*Math.cos(theta2-phi)*Math.sin(phi) + v1*Math.sin(theta1-phi)*Math.sin(phi+0.5*Math.PI)
                b2[2] = v1*Math.cos(theta1-phi)*Math.cos(phi) + v2*Math.sin(theta2-phi)*Math.cos(phi+0.5*Math.PI)
                b2[3] = v1*Math.cos(theta1-phi)*Math.sin(phi) + v2*Math.sin(theta2-phi)*Math.sin(phi+0.5*Math.PI)
                // evolve trajectory post-collision (for dt2)
                const dt2 = dt - dt1
                // ball 1
                b1[0] = b1[0] + b1[2]*dt2 + 0.5*ax*dt2*dt2 // x
                b1[1] = b1[1] + b1[3]*dt2 + 0.5*ay*dt2*dt2 // y
                b1[2] = b1[2] + ax*dt2 // vx
                b1[3] = b1[3] + ay*dt2 // vy
                // ball 2
                b2[0] = b2[0] + b2[2]*dt2 + 0.5*ax*dt2*dt2 // x
                b2[1] = b2[1] + b2[3]*dt2 + 0.5*ay*dt2*dt2 // y
                b2[2] = b2[2] + ax*dt2 // vx
                b2[3] = b2[3] + ay*dt2 // vy
            }
        }
    }
    
    updateBallStates(dt) {
        this.updateBallisticTrajectory(dt)
        this.deleteOffScreenBalls()
        this.handleBallBallCollisions(dt)
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