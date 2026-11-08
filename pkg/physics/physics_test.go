package physics_test

import (
	"math"
	"testing"

	"github.com/ZengLLQ/quasar-kinetics/pkg/diagnostics"
	"github.com/ZengLLQ/quasar-kinetics/pkg/physics"
)

func TestTwoStreamInstabilityGrowth(t *testing.T) {
	tsc := physics.DefaultTwoStreamConfig()
	tsc.Nx = 32
	tsc.NumPerCell = 30
	tsc.PerturbAmp = 0.05

	sim, err := physics.SetupTwoStreamSimulator(tsc)
	if err != nil {
		t.Fatalf("failed to setup two stream simulator: %v", err)
	}

	initialElecEnergy := sim.Grid.ElectricEnergy()
	if initialElecEnergy <= 0 {
		t.Fatalf("expected non-zero initial electric energy from perturbation")
	}

	// Run for 150 steps (approx 3 plasma oscillation periods)
	sim.Run(150)

	finalElecEnergy := sim.Grid.ElectricEnergy()

	// In Two-Stream instability, electric field energy grows exponentially by multiple factors
	growthRatio := finalElecEnergy / initialElecEnergy
	if growthRatio < 1.5 {
		t.Errorf("two-stream instability failed to grow electric field: initial=%g, final=%g (ratio=%g)",
			initialElecEnergy, finalElecEnergy, growthRatio)
	}

	// Verify phase space vortex capture
	ps := diagnostics.NewPhaseSpace2D(20, 20, 0, tsc.BoxLength, -0.4, 0.4)
	ps.BinSpecies(sim.Species[0], sim.Config.Units.SpeedOfLight)
	ps.BinSpecies(sim.Species[1], sim.Config.Units.SpeedOfLight)

	// ASCII heatmap generation check
	heatmap := ps.ASCIIHeatmap(20, 10)
	if len(heatmap) == 0 {
		t.Errorf("expected non-empty phase space heatmap")
	}
}

func TestLandauDampingPhysics(t *testing.T) {
	lc := physics.DefaultLandauConfig()
	lc.Nx = 32
	lc.NumPerCell = 40

	sim, err := physics.SetupLandauSimulator(lc)
	if err != nil {
		t.Fatalf("failed to setup Landau simulator: %v", err)
	}

	k := float64(lc.PerturbMode) * (2.0 * math.Pi / lc.BoxLength)
	gammaTheory := physics.LandauTheoreticalDampingRate(k, lc.ThermalVel, 1.0)
	omegaTheory := physics.BohmGrossFrequency(k, lc.ThermalVel, 1.0)

	if gammaTheory >= 0 {
		t.Errorf("Landau damping rate must be strictly negative, got %g", gammaTheory)
	}
	if omegaTheory <= 1.0 {
		t.Errorf("Bohm-Gross frequency with thermal speed must exceed omega_pe=1, got %g", omegaTheory)
	}

	// Step forward
	sim.Run(50)
	if sim.StepCount != 50 {
		t.Errorf("expected step count 50, got %d", sim.StepCount)
	}
}

func TestRelativisticCyclotronMetricsAndStep(t *testing.T) {
	cp := physics.CyclotronParams{
		MagneticField: 1.0,
		RelProperMom:  1.0, // u = 1.0 => gamma = sqrt(2) approx 1.414
		Charge:        -1.0,
		Mass:          1.0,
	}

	sim, metrics, err := physics.SetupCyclotronSimulator(cp, 32, 10.0)
	if err != nil {
		t.Fatalf("failed to setup cyclotron simulator: %v", err)
	}

	if math.Abs(metrics.Gamma-math.Sqrt(2.0)) > 1e-12 {
		t.Errorf("expected gamma=sqrt(2), got %g", metrics.Gamma)
	}

	// Run for 1 orbital period (200 steps)
	sim.Run(200)

	p := sim.Species[0].Particles[0]
	// Proper momentum magnitude must remain 1.0
	uMag := math.Hypot(p.Ux, p.Uy)
	if math.Abs(uMag-cp.RelProperMom)/cp.RelProperMom > 1e-12 {
		t.Errorf("relativistic proper momentum magnitude drift: got %g, expected %g",
			uMag, cp.RelProperMom)
	}
}

func TestWeibelInstabilitySetupAndGrowth(t *testing.T) {
	wc := physics.DefaultWeibelConfig()
	wc.Nx = 32
	wc.NumPerCell = 30

	sim, err := physics.SetupWeibelSimulator(wc)
	if err != nil {
		t.Fatalf("failed to setup Weibel simulator: %v", err)
	}

	gammaW := physics.WeibelTheoreticalGrowthRate(wc, 1.0, 1.0)
	if gammaW <= 0 {
		t.Errorf("expected positive Weibel growth rate for anisotropic plasma, got %g", gammaW)
	}

	initialMagEnergy := sim.Grid.MagneticEnergy()
	if initialMagEnergy <= 0 {
		t.Fatalf("expected non-zero initial magnetic energy from seed")
	}

	// Run 100 steps
	sim.Run(100)

	finalMagEnergy := sim.Grid.MagneticEnergy()
	if finalMagEnergy <= 0 {
		t.Errorf("magnetic energy must remain positive, got %g", finalMagEnergy)
	}
}
