package physics

import (
	"math"

	"github.com/ZengLLQ/quasar-kinetics/pkg/constants"
	"github.com/ZengLLQ/quasar-kinetics/pkg/engine"
	"github.com/ZengLLQ/quasar-kinetics/pkg/particles"
)

// WeibelConfig specifies parameters for the electromagnetic Weibel filamentation instability.
type WeibelConfig struct {
	Nx          int     // Grid cells
	NumPerCell  int     // Particles per cell
	VthX        float64 // Thermal spread along x (low)
	VthY        float64 // Thermal spread along y (high, anisotropic)
	BoxLength   float64 // Domain length
	PerturbAmpB float64 // Small initial magnetic field seed
}

// DefaultWeibelConfig returns standard parameters for magnetic filamentation growth.
func DefaultWeibelConfig() WeibelConfig {
	return WeibelConfig{
		Nx:          64,
		NumPerCell:  50,
		VthX:        0.05,
		VthY:        0.30, // Strongly anisotropic Ty >> Tx
		BoxLength:   2.0 * math.Pi,
		PerturbAmpB: 1e-4,
	}
}

// WeibelTheoreticalGrowthRate calculates the linear growth rate of transverse magnetic field:
//
//	gamma_W approx omega_pe * (vth_y / c) * sqrt( (Ty - Tx) / Ty )
func WeibelTheoreticalGrowthRate(wc WeibelConfig, omegaPe, c float64) float64 {
	tx := wc.VthX * wc.VthX
	ty := wc.VthY * wc.VthY
	if ty <= tx {
		return 0.0
	}
	anisotropy := (ty - tx) / ty
	return omegaPe * (wc.VthY / c) * math.Sqrt(anisotropy)
}

// SetupWeibelSimulator initializes a Simulator for the Weibel instability.
func SetupWeibelSimulator(wc WeibelConfig) (*engine.Simulator, error) {
	cfg := engine.DefaultConfig(wc.Nx, wc.BoxLength)
	cfg.Dt = 0.02

	sim, err := engine.NewSimulator(cfg)
	if err != nil {
		return nil, err
	}

	totalParticles := wc.Nx * wc.NumPerCell
	macroWeight := (wc.BoxLength / float64(totalParticles)) / cfg.Dx

	electrons := particles.NewSpecies(0, "anisotropic_electrons", -1.0, 1.0, totalParticles)
	ions := particles.NewSpecies(1, "ions", 1.0, 1836.0, totalParticles)

	for i := 0; i < totalParticles; i++ {
		x := float64(i) * (wc.BoxLength / float64(totalParticles))

		// Box-Muller anisotropic velocity generation
		u1 := (float64(i%127) + 0.5) / 127.0
		u2 := (float64((i*37)%127) + 0.5) / 127.0
		r1 := math.Sqrt(-2.0 * math.Log(u1))
		theta1 := 2.0 * math.Pi * u2
		vx := wc.VthX * r1 * math.Cos(theta1)
		vy := wc.VthY * r1 * math.Sin(theta1)

		gamma := constants.RelativisticGamma(vx, vy, 0, cfg.Units.SpeedOfLight)

		electrons.AddParticle(particles.Particle{
			X:      x,
			Ux:     gamma * vx,
			Uy:     gamma * vy,
			Weight: macroWeight,
		})

		ions.AddParticle(particles.Particle{
			X:      x,
			Ux:     0.0,
			Weight: macroWeight,
		})
	}

	// Seed small magnetic fluctuation in Bz to trigger mode 1
	k1 := 2.0 * math.Pi / wc.BoxLength
	for i := 0; i < sim.Grid.Nx; i++ {
		x := float64(i) * sim.Grid.Dx
		sim.Grid.Bz[sim.Grid.Nghost+i] = wc.PerturbAmpB * math.Sin(k1*x)
	}
	sim.Grid.ApplyBoundaries()

	sim.AddSpecies(electrons)
	sim.AddSpecies(ions)

	sim.InitializeElectrostatic()
	return sim, nil
}
