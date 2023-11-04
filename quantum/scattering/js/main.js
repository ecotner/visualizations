import { BallGame } from "./balls.js"

// init game
let game = new BallGame(300, 300, 20)

// attach game screen to DOM
let canvas = game.cm.canvas
canvas.style.border = "solid black 1px"
let body = document.getElementById("body")
body.append(canvas)

// start game loop
game.run()
