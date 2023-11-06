import { BallGame } from "./balls.js"

// init game
const game = new BallGame(600, 600, 20)

// attach game screen to DOM
const canvas = game.cm.canvas
canvas.style.border = "solid black 1px"
const body = document.getElementById("body")
body.append(canvas)

// set up dropdown for switching between ball/wall creation
const dropdown = document.createElement("select")
for (let x of ["balls", "walls"]) {
    const option = document.createElement("option")
    option.innerHTML = `Create ${x}`
    option.value = x
    dropdown.options.add(option)
}
dropdown.addEventListener("change", (event) => {
    game.objCreationState = event.target.value
})
body.append(dropdown)

// start game loop
game.run()
