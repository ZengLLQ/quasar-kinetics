package diagnostics

import (
	"math"

	"github.com/ZengLLQ/quasar-kinetics/pkg/fields"
)

// SpectralAnalyzer performs spatial Fourier decomposition of electromagnetic field profiles.
type SpectralAnalyzer struct{}

// NewSpectralAnalyzer creates an initialized SpectralAnalyzer.
func NewSpectralAnalyzer() *SpectralAnalyzer {
	return &SpectralAnalyzer{}
}

// ModeAmplitude represents an individual spatial Fourier mode.
type ModeAmplitude struct {
	ModeNumber int     `json:"mode_number"`
	Wavenumber float64 `json:"wavenumber"`
	Real       float64 `json:"real"`
	Imag       float64 `json:"imag"`
	Amplitude  float64 `json:"amplitude"`
}

// AnalyzeExModes calculates spatial Fourier amplitudes of longitudinal electric field Ex.
func (sa *SpectralAnalyzer) AnalyzeExModes(g *fields.Grid1D) []ModeAmplitude {
	nx := g.Nx
	ng := g.Nghost
	numModes := nx / 2

	modes := make([]ModeAmplitude, numModes)
	dk := (2.0 * math.Pi) / g.Lx

	for m := 0; m < numModes; m++ {
		k := float64(m) * dk
		var realSum, imagSum float64

		for j := 0; j < nx; j++ {
			x := float64(j) * g.Dx
			theta := k * x
			val := g.Ex[ng+j]
			realSum += val * math.Cos(theta)
			imagSum -= val * math.Sin(theta)
		}

		norm := 2.0 / float64(nx)
		if m == 0 {
			norm = 1.0 / float64(nx)
		}
		realSum *= norm
		imagSum *= norm

		amp := math.Hypot(realSum, imagSum)
		modes[m] = ModeAmplitude{
			ModeNumber: m,
			Wavenumber: k,
			Real:       realSum,
			Imag:       imagSum,
			Amplitude:  amp,
		}
	}

	return modes
}

// DominantWaveMode finds the non-zero Fourier mode with largest electric field amplitude.
func (sa *SpectralAnalyzer) DominantWaveMode(g *fields.Grid1D) ModeAmplitude {
	modes := sa.AnalyzeExModes(g)
	if len(modes) <= 1 {
		return ModeAmplitude{}
	}

	maxAmp := -1.0
	var dom ModeAmplitude
	for m := 1; m < len(modes); m++ {
		if modes[m].Amplitude > maxAmp {
			maxAmp = modes[m].Amplitude
			dom = modes[m]
		}
	}
	return dom
}
