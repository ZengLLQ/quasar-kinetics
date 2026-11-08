package deposition

import "math"

// ShapeFunctionOrder defines the interpolation order for particle-mesh operations.
type ShapeFunctionOrder int

const (
	// ShapeNGP represents 0th-order Nearest-Grid-Point weighting.
	ShapeNGP ShapeFunctionOrder = 0
	// ShapeCIC represents 1st-order Cloud-In-Cell (linear B-spline) weighting.
	ShapeCIC ShapeFunctionOrder = 1
	// ShapeTSC represents 2nd-order Triangular-Shaped-Cloud (quadratic B-spline) weighting.
	ShapeTSC ShapeFunctionOrder = 2
)

// WeightsCIC computes the 1st-order linear Cloud-in-Cell weights onto adjacent nodes.
// u is the continuous coordinate in units of grid spacing: u = (x - x0) / dx.
// It returns base node index i = floor(u), weight w0 on node i, and w1 on node i+1.
func WeightsCIC(u float64) (i int, w0, w1 float64) {
	i = int(math.Floor(u))
	frac := u - float64(i)
	w0 = 1.0 - frac
	w1 = frac
	return i, w0, w1
}

// WeightsTSC computes the 2nd-order quadratic B-spline weights onto three adjacent nodes:
// i-1, i, and i+1.
func WeightsTSC(u float64) (i int, wm1, w0, wp1 float64) {
	i = int(math.Floor(u + 0.5))
	delta := u - float64(i) // delta in [-0.5, 0.5]

	wm1 = 0.5 * (0.5 - delta) * (0.5 - delta)
	w0 = 0.75 - delta*delta
	wp1 = 0.5 * (0.5 + delta) * (0.5 + delta)
	return i, wm1, w0, wp1
}
