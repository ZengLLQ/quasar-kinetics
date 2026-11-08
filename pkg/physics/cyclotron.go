package physics

import (
	"math"

	"github.com/ZengLLQ/quasar-kinetics/pkg/engine"
	"github.com/ZengLLQ/quasar-kinetics/pkg/particles"
)

// CyclotronParams defines settings for relativistic cyclotron motion.
type CyclotronParams struct {
	MagneticField float64 // Guide field Bz [T or normalized]
	RelProperMom  float64 // Relativistic proper velocity u0 = gamma * v0
	Charge        float64 // Particle charge
	Mass          float64 // Rest mass
}

// CyclotronAnalyticalMetrics calculates exact relativistic gyration metrics.
type CyclotronAnalyticalMetrics struct {
	Gamma        float64 // Lorentz factor
	Velocity     float64 // Physical velocity v0
	OmegaC       float64 // Relativistic gyrofrequency omega_c = q * B / (gamma * m)
	LarmorRadius float64 // Gyroradius r_L = gamma * m * v0 / (q * B)
	Period       float64 // Orbital period T = 2 * pi / omega_c
}

// ComputeCyclotronMetrics evaluates theoretical relativistic cyclotron values.
func ComputeCyclotronMetrics(cp CyclotronParams, c float64) CyclotronAnalyticalMetrics {
	gamma := math.Sqrt(1.0 + (cp.RelProperMom*cp.RelProperMom)/(c*c))
	v0 := cp.RelProperMom / gamma

	omegaC := math.Abs(cp.Charge*cp.MagneticField) / (gamma * cp.Mass)
	rL := (gamma * cp.Mass * v0) / math.Abs(cp.Charge*cp.MagneticField)
	period := 2.0 * math.Pi / omegaC

	return CyclotronAnalyticalMetrics{
		Gamma:        gamma,
		Velocity:     v0,
		OmegaC:       omegaC,
		LarmorRadius: rL,
		Period:       period,
	}
}

// SetupCyclotronSimulator initializes a simulation tracking a relativistic particle in a magnetic field.
func SetupCyclotronSimulator(cp CyclotronParams, nx int, lx float64) (*engine.Simulator, CyclotronAnalyticalMetrics, error) {
	cfg := engine.DefaultConfig(nx, lx)
	metrics := ComputeCyclotronMetrics(cp, cfg.Units.SpeedOfLight)
	cfg.Dt = metrics.Period / 200.0 // 200 steps per orbital period

	sim, err := engine.NewSimulator(cfg)
	if err != nil {
		return nil, metrics, err
	}

	// Apply uniform magnetic guide field Bz across all grid cells
	for i := 0; i < sim.Grid.TotalX; i++ {
		sim.Grid.Bz[i] = cp.MagneticField
	}
	sim.Grid.ApplyBoundaries()

	// Single relativistic test particle initialized at top of orbit (weight = 0 for pure test particle)
	sp := particles.NewSpecies(0, "rel_electron", cp.Charge, cp.Mass, 1)
	sp.AddParticle(particles.Particle{
		X:      0.5 * lx,
		Y:      metrics.LarmorRadius,
		Ux:     cp.RelProperMom,
		Uy:     0.0,
		Uz:     0.0,
		Weight: 0.0,
	})

	sim.AddSpecies(sp)
	return sim, metrics, nil
}
