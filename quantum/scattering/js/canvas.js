export class CanvasManager {
    /**
     * Allows for high-level manipulation of the canvas
     * @param {number} height The height of the canvas element in px
     * @param {number} width The width of the canvas element in px
     */
    constructor(height, width) {
        // set up canvas element
        this.canvas = document.createElement("canvas")
        this.canvas.id = "canvas"
        this.height = height
        this.width = width
        this.ctx = this.canvas.getContext("2d")
        // transform canvas coordinates so x increases left -> right and y increases
        // bottom -> top, starting with (0, 0) in the bottom left corner
        this.ctx.transform(1, 0, 0, -1, 0, height)
        this.ctx.save() // save initial blank state if we need to reset
        // create a buffer to use for writing image data into canvas; pre-set opacity to 100%
        this.imgDataBuffer = this.ctx.createImageData(this.width, this.height)
        for (let i = 0; i < height * width; i++) {
            this.imgDataBuffer.data[4*i + 3] = 255
        }
    }

    get height() {
        return this.canvas.height
    }
    
    get width() {
        return this.canvas.width
    }

    set height(val) {
        this.canvas.height = val
    }

    set width(val) {
        this.canvas.width = val
    }

    /**
     * renders an array of floats as a rasterized monochrome image
     * @param {Float32Array} array An array of float values, scaled to be between 0 and 1
     */
    render_raster_float(array) {
        const data = this.imgDataBuffer.data
        for (let i = 0; i < array.length; i ++) {
            data[4*i] = (255 * array[i]) | 0 // monochrome red
        }
        this.ctx.putImageData(this.imgDataBuffer, 0, 0)
    }

    /**
     * renders an array of u8 integers as a rasterized monochrome image
     * @param {Uint8Array} array An array of u8 values
     */
    render_raster_u8(array) {
        const data = this.imgDataBuffer.data
        for (let i = 0; i < array.length; i ++) {
            data[4*i] = array[i] // monochrome red
        }
        this.ctx.putImageData(this.imgDataBuffer, 0, 0)
    }
}