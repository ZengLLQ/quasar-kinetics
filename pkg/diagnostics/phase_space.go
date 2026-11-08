package diagnostics

import (
	"fmt"
	"math"
	"strings"

	"github.com/ZengLLQ/quasar-kinetics/pkg/particles"
)

// PhaseSpace2D represents a discretized 2D histogram of phase space distribution f(x, p_x).
type PhaseSpace2D struct {
	Nx     int         // Spatial bin count
	Nv     int         // Velocity/momentum bin count
	XMin   float64     // Minimum x
	XMax   float64     // Maximum x
	VMin   float64     // Minimum v/u
	VMax   float64     // Maximum v/u
	Dx     float64     // Spatial bin width
	Dv     float64     // Velocity bin width
	Counts [][]float64 // [nx][nv] weighted particle density
}

// NewPhaseSpace2D creates an empty phase space histogram.
func NewPhaseSpace2D(nx, nv int, xMin, xMax, vMin, vMax float64) *PhaseSpace2D {
	dx := (xMax - xMin) / float64(nx)
	dv := (vMax - vMin) / float64(nv)

	counts := make([][]float64, nx)
	for i := 0; i < nx; i++ {
		counts[i] = make([]float64, nv)
	}

	return &PhaseSpace2D{
		Nx:     nx,
		Nv:     nv,
		XMin:   xMin,
		XMax:   xMax,
		VMin:   vMin,
		VMax:   vMax,
		Dx:     dx,
		Dv:     dv,
		Counts: counts,
	}
}

// BinSpecies deposits particle weights of a species into the phase space histogram.
func (ps *PhaseSpace2D) BinSpecies(sp *particles.Species, c float64) {
	for i := range sp.Particles {
		p := &sp.Particles[i]
		vx, _, _ := p.Velocity(c)

		// Spatial bin
		uX := (p.X - ps.XMin) / ps.Dx
		xIdx := int(math.Floor(uX))
		if xIdx < 0 || xIdx >= ps.Nx {
			continue
		}

		// Velocity bin
		uV := (vx - ps.VMin) / ps.Dv
		vIdx := int(math.Floor(uV))
		if vIdx < 0 || vIdx >= ps.Nv {
			continue
		}

		ps.Counts[xIdx][vIdx] += p.Weight
	}
}

// ASCIIHeatmap renders an ASCII art visualization of phase space distribution f(x, v_x).
func (ps *PhaseSpace2D) ASCIIHeatmap(width, height int) string {
	if width <= 0 || width > ps.Nx {
		width = ps.Nx
	}
	if height <= 0 || height > ps.Nv {
		height = ps.Nv
	}

	chars := []rune(" .:-=+*#%@")
	numChars := len(chars)

	// Find max density
	var maxVal float64
	for i := 0; i < ps.Nx; i++ {
		for j := 0; j < ps.Nv; j++ {
			if ps.Counts[i][j] > maxVal {
				maxVal = ps.Counts[i][j]
			}
		}
	}
	if maxVal == 0 {
		maxVal = 1.0
	}

	var sb strings.Builder
	sb.WriteString(fmt.Sprintf("Phase Space f(x, vx) [x: [%.2f, %.2f], vx: [%.2f, %.2f]]\n",
		ps.XMin, ps.XMax, ps.VMin, ps.VMax))
	sb.WriteString(strings.Repeat("-", width+12) + "\n")

	// Render from high velocity (top) to low velocity (bottom)
	for row := height - 1; row >= 0; row-- {
		vVal := ps.VMin + (float64(row)+0.5)*(ps.VMax-ps.VMin)/float64(height)
		sb.WriteString(fmt.Sprintf("%+6.2f | ", vVal))

		for col := 0; col < width; col++ {
			xBin := col * ps.Nx / width
			vBin := row * ps.Nv / height

			val := ps.Counts[xBin][vBin]
			charIdx := int((val / maxVal) * float64(numChars-1))
			if charIdx >= numChars {
				charIdx = numChars - 1
			} else if charIdx < 0 {
				charIdx = 0
			}
			sb.WriteRune(chars[charIdx])
		}
		sb.WriteString(" |\n")
	}
	sb.WriteString(strings.Repeat("-", width+12) + "\n")
	sb.WriteString("       | " + strings.Repeat("^", width) + " |\n")
	sb.WriteString(fmt.Sprintf("       x: %.2f%s%.2f\n",
		ps.XMin, strings.Repeat(" ", width-12), ps.XMax))

	return sb.String()
}
