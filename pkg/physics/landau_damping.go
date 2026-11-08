package physics

import (
	"math"

	"github.com/ZengLLQ/quasar-kinetics/pkg/constants"
	"github.com/ZengLLQ/quasar-kinetics/pkg/engine"
	"github.com/ZengLLQ/quasar-kinetics/pkg/particles"
)

// LandauConfig holds parameters for a Landau damping simulation.
type LandauConfig struct {
	Nx          int     // Grid cells
	NumPerCell  int     // Particles per cell
	ThermalVel  float64 // Electron thermal speed v_th = sqrt(k_B * T / m)
	PerturbAmp  float64 // Perturbation amplitude alpha
	PerturbMode int     // Perturbed mode index m
	BoxLength   float64 // Domain length Lx
}

// DefaultLandauConfig returns standard parameters for collisionless Landau damping.
func DefaultLandauConfig() LandauConfig {
	return LandauConfig{
		Nx:          64,
		NumPerCell:  100,
		ThermalVel:  0.2,
		PerturbAmp:  0.05,
		PerturbMode: 1,
		BoxLength:   2.0 * math.Pi,
	}
}

// LandauTheoreticalDampingRate calculates linear Landau damping rate:
//
//	gamma_L = - sqrt(pi / 8) * (omega_pe / (k * lambda_D)^3) * exp( - 1 / (2 * (k * lambda_D)^2) - 1.5 )
func LandauTheoreticalDampingRate(k, vTh, omegaPe float64) float64 {
	lambdaD := vTh / omegaPe
	kLambdaD := k * lambdaD
	if kLambdaD <= 0 {
		return 0
	}

	coeff := -math.Sqrt(math.Pi/8.0) * omegaPe / math.Pow(kLambdaD, 3)
	arg := -1.0/(2.0*kLambdaD*kLambdaD) - 1.5
	return coeff * math.Exp(arg)
}

// BohmGrossFrequency computes the Langmuir wave frequency with thermal dispersion:
//
//	omega^2 = omega_pe^2 + 3 * k^2 * v_th^2
func BohmGrossFrequency(k, vTh, omegaPe float64) float64 {
	return math.Sqrt(omegaPe*omegaPe + 3.0*k*k*vTh*vTh)
}

// SetupLandauSimulator initializes a Simulator for observing Landau damping of an electrostatic wave.
func SetupLandauSimulator(lc LandauConfig) (*engine.Simulator, error) {
	cfg := engine.DefaultConfig(lc.Nx, lc.BoxLength)
	cfg.Dt = 0.02

	sim, err := engine.NewSimulator(cfg)
	if err != nil {
		return nil, err
	}

	totalParticles := lc.Nx * lc.NumPerCell
	macroWeight := (lc.BoxLength / float64(totalParticles)) / cfg.Dx
	kPerturb := float64(lc.PerturbMode) * (2.0 * math.Pi / lc.BoxLength)

	electrons := particles.NewSpecies(0, "electrons", -1.0, 1.0, totalParticles)
	ions := particles.NewSpecies(1, "ions", 1.0, 1836.0, totalParticles)

	for i := 0; i < totalParticles; i++ {
		xBase := float64(i) * (lc.BoxLength / float64(totalParticles))

		// Sinusoidal density perturbation: x' = x + (alpha / k) * sin(k * x)
		xPerturbed := xBase + (lc.PerturbAmp/kPerturb)*math.Sin(kPerturb*xBase)
		xWrapped := math.Mod(xPerturbed, lc.BoxLength)
		if xWrapped < 0 {
			xWrapped += lc.BoxLength
		}

		// Box-Muller Gaussian velocity sampling
		u1 := (float64(i%137) + 0.5) / 137.0
		u2 := (float64((i*43)%137) + 0.5) / 137.0
		r := math.Sqrt(-2.0 * math.Log(u1))
		theta := 2.0 * math.Pi * u2
		vx := lc.ThermalVel * r * math.Cos(theta)
		vy := lc.ThermalVel * r * math.Sin(theta)

		gamma := constants.RelativisticGamma(vx, vy, 0, cfg.Units.SpeedOfLight)

		electrons.AddParticle(particles.Particle{
			X:      xWrapped,
			Ux:     gamma * vx,
			Uy:     gamma * vy,
			Weight: macroWeight,
		})

		// Uniform neutralizing background ions
		ions.AddParticle(particles.Particle{
			X:      xBase,
			Ux:     0.0,
			Weight: macroWeight,
		})
	}

	sim.AddSpecies(electrons)
	sim.AddSpecies(ions)

	sim.InitializeElectrostatic()
	return sim, nil
}
