package physics

import (
	"math"

	"github.com/ZengLLQ/quasar-kinetics/pkg/constants"
	"github.com/ZengLLQ/quasar-kinetics/pkg/engine"
	"github.com/ZengLLQ/quasar-kinetics/pkg/particles"
)

// TwoStreamConfig defines parameters for a two-stream instability simulation.
type TwoStreamConfig struct {
	Nx          int     // Grid cells
	NumPerCell  int     // Particles per cell per beam
	DriftVel    float64 // Beam drift velocity v0 (e.g. 0.2)
	ThermalVel  float64 // Beam thermal spread v_th (e.g. 0.01)
	PerturbAmp  float64 // Initial spatial perturbation amplitude (e.g. 0.01)
	PerturbMode int     // Wavenumber mode index to perturb (e.g. 1)
	BoxLength   float64 // Total physical domain length Lx
}

// DefaultTwoStreamConfig returns standard benchmark settings.
func DefaultTwoStreamConfig() TwoStreamConfig {
	// For v0 = 0.2 and omega_pe = 1, unstable modes satisfy k * v0 < omega_pe
	// Choosing Lx = 2 * pi allows mode 1 (k = 1) to have k * v0 = 0.2 < 1.0 (unstable!)
	return TwoStreamConfig{
		Nx:          64,
		NumPerCell:  20,
		DriftVel:    0.2,
		ThermalVel:  0.005,
		PerturbAmp:  0.02,
		PerturbMode: 1,
		BoxLength:   2.0 * math.Pi,
	}
}

// TwoStreamAnalyticalGrowthRate computes the maximum linear growth rate gamma from Vlasov cold beam theory:
//
//	gamma_max = omega_pe / (2 * sqrt(2)) approx 0.3535 * omega_pe
func TwoStreamAnalyticalGrowthRate(omegaPe float64) float64 {
	return omegaPe / (2.0 * math.Sqrt(2.0))
}

// SetupTwoStreamSimulator initializes a Simulator pre-configured with two counter-streaming electron beams
// and a neutralizing background ion population.
func SetupTwoStreamSimulator(tsc TwoStreamConfig) (*engine.Simulator, error) {
	cfg := engine.DefaultConfig(tsc.Nx, tsc.BoxLength)
	cfg.Dt = 0.02

	sim, err := engine.NewSimulator(cfg)
	if err != nil {
		return nil, err
	}

	totalParticlesPerBeam := tsc.Nx * tsc.NumPerCell
	macroWeight := (tsc.BoxLength / float64(totalParticlesPerBeam)) / cfg.Dx

	kPerturb := float64(tsc.PerturbMode) * (2.0 * math.Pi / tsc.BoxLength)

	// Beam 1: drifting with +v0
	beam1 := particles.NewSpecies(0, "beam_plus", -1.0, 1.0, totalParticlesPerBeam)
	// Beam 2: drifting with -v0
	beam2 := particles.NewSpecies(1, "beam_minus", -1.0, 1.0, totalParticlesPerBeam)
	// Neutralizing ions: stationary
	ions := particles.NewSpecies(2, "ions", 1.0, 1836.0, totalParticlesPerBeam*2)

	for i := 0; i < totalParticlesPerBeam; i++ {
		// Uniform spatial placement
		xBase := float64(i) * (tsc.BoxLength / float64(totalParticlesPerBeam))

		// Apply sinusoidal spatial perturbation to seed mode
		xPerturbed := xBase + (tsc.PerturbAmp/kPerturb)*math.Sin(kPerturb*xBase)
		xWrapped := math.Mod(xPerturbed, tsc.BoxLength)
		if xWrapped < 0 {
			xWrapped += tsc.BoxLength
		}

		// Gaussian thermal dispersion (Box-Muller)
		u1 := (float64(i%97) + 0.5) / 97.0
		u2 := (float64((i*31)%97) + 0.5) / 97.0
		r := math.Sqrt(-2.0 * math.Log(u1))
		theta := 2.0 * math.Pi * u2
		vTh1 := tsc.ThermalVel * r * math.Cos(theta)
		vTh2 := tsc.ThermalVel * r * math.Sin(theta)

		gamma1 := constants.RelativisticGamma(tsc.DriftVel+vTh1, 0, 0, cfg.Units.SpeedOfLight)
		gamma2 := constants.RelativisticGamma(-tsc.DriftVel+vTh2, 0, 0, cfg.Units.SpeedOfLight)

		beam1.AddParticle(particles.Particle{
			X:      xWrapped,
			Ux:     gamma1 * (tsc.DriftVel + vTh1),
			Weight: macroWeight,
		})

		beam2.AddParticle(particles.Particle{
			X:      xWrapped,
			Ux:     gamma2 * (-tsc.DriftVel + vTh2),
			Weight: macroWeight,
		})

		// 2 neutralizing ions per electron pair
		ions.AddParticle(particles.Particle{
			X:      xBase,
			Ux:     0.0,
			Weight: macroWeight,
		})
		ions.AddParticle(particles.Particle{
			X:      xBase + 0.5*(tsc.BoxLength/float64(totalParticlesPerBeam)),
			Ux:     0.0,
			Weight: macroWeight,
		})
	}

	sim.AddSpecies(beam1)
	sim.AddSpecies(beam2)
	sim.AddSpecies(ions)

	sim.InitializeElectrostatic()
	return sim, nil
}
