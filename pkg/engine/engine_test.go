package engine_test

import (
	"math"
	"testing"

	"github.com/ZengLLQ/quasar-kinetics/pkg/engine"
	"github.com/ZengLLQ/quasar-kinetics/pkg/particles"
)

func TestConfigValidation(t *testing.T) {
	cfg := engine.DefaultConfig(64, 10.0)
	if err := cfg.Validate(); err != nil {
		t.Fatalf("default config should be valid, got: %v", err)
	}

	// Invalid Nx
	cfgBadNx := *cfg
	cfgBadNx.Nx = 1
	if err := cfgBadNx.Validate(); err == nil {
		t.Errorf("expected error for Nx <= 2")
	}

	// Invalid CFL
	cfgBadCFL := *cfg
	cfgBadCFL.CFL = 1.5
	if err := cfgBadCFL.Validate(); err == nil {
		t.Errorf("expected error for CFL > 1.0")
	}
}

func TestTilePartitionBinning(t *testing.T) {
	lx := 10.0
	x0 := 0.0
	numTiles := 5
	tiler := engine.NewTilePartition(numTiles, x0, lx)

	sp := particles.NewSpecies(0, "test", -1.0, 1.0, 10)
	// Place particles at 1.0, 3.0, 5.0, 7.0, 9.0
	// With width = 2.0:
	// 1.0 in tile 0 [0, 2)
	// 3.0 in tile 1 [2, 4)
	// 5.0 in tile 2 [4, 6)
	// 7.0 in tile 3 [6, 8)
	// 9.0 in tile 4 [8, 10)
	positions := []float64{1.0, 3.0, 5.0, 7.0, 9.0}
	for _, p := range positions {
		sp.AddParticle(particles.Particle{X: p, Weight: 1.0})
	}

	tiler.BinSpecies(sp)

	for tileIdx := 0; tileIdx < numTiles; tileIdx++ {
		if len(tiler.Tiles[tileIdx].Indices) != 1 {
			t.Errorf("tile %d expected 1 particle, got %d",
				tileIdx, len(tiler.Tiles[tileIdx].Indices))
		}
	}
}

func TestSimulatorEnergyConservationDrift(t *testing.T) {
	cfg := engine.DefaultConfig(64, 2.0*math.Pi)
	cfg.Dt = 0.02
	cfg.TotalSteps = 100
	cfg.DiagnosticInterval = 10
	cfg.NumWorkers = 2

	sim, err := engine.NewSimulator(cfg)
	if err != nil {
		t.Fatalf("failed to create simulator: %v", err)
	}

	// Add two counter-streaming electron beams so net current is zero: <J_x> = 0
	electrons := particles.NewSpecies(0, "electrons", -1.0, 1.0, 200)
	dx := cfg.Dx
	nParticles := 200
	for i := 0; i < nParticles; i++ {
		x := float64(i) * (sim.Grid.Lx / float64(nParticles))
		ux := 0.1
		if i%2 == 1 {
			ux = -0.1
		}
		electrons.AddParticle(particles.Particle{
			X:      x,
			Ux:     ux,
			Weight: (sim.Grid.Lx / float64(nParticles)) / dx,
		})
	}
	sim.AddSpecies(electrons)

	// Add uniform background ions for charge neutrality
	ions := particles.NewSpecies(1, "ions", 1.0, 1836.0, 200)
	for i := 0; i < nParticles; i++ {
		x := float64(i) * (sim.Grid.Lx / float64(nParticles))
		ions.AddParticle(particles.Particle{
			X:      x,
			Ux:     0.0,
			Weight: (sim.Grid.Lx / float64(nParticles)) / dx,
		})
	}
	sim.AddSpecies(ions)

	sim.InitializeElectrostatic()
	eKin0, eElec0, eMag0, eTot0 := sim.TotalEnergy()

	if eTot0 <= 0 {
		t.Fatalf("initial total energy should be positive, got: %g", eTot0)
	}
	_ = eKin0
	_ = eElec0
	_ = eMag0

	// Run 100 steps
	sim.Run(100)

	_, _, _, eTotFinal := sim.TotalEnergy()
	drift := math.Abs(eTotFinal-eTot0) / eTot0

	// Total energy conservation should be within 2.5% for 100 leapfrog PIC steps
	if drift > 0.025 {
		t.Errorf("total energy drift exceeded threshold: drift=%g (initial=%g, final=%g)",
			drift, eTot0, eTotFinal)
	}

	if len(sim.History) < 5 {
		t.Errorf("expected at least 5 diagnostic records, got %d", len(sim.History))
	}
}
