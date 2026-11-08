package particles_test

import (
	"math"
	"testing"

	"github.com/ZengLLQ/quasar-kinetics/pkg/constants"
	"github.com/ZengLLQ/quasar-kinetics/pkg/fields"
	"github.com/ZengLLQ/quasar-kinetics/pkg/particles"
)

func TestBorisPureMagneticConservation(t *testing.T) {
	// In pure magnetic field (E=0), kinetic energy and proper momentum magnitude must be exactly conserved
	units := constants.NewNormalizedUnits(1.0, 1.0)
	pusher := particles.NewBorisPusher(units)

	p := particles.Particle{
		X:  0.0,
		Y:  0.0,
		Z:  0.0,
		Ux: 1.5, // Highly relativistic (u/c = 1.5 => gamma = sqrt(1 + 2.25) = 1.80277)
		Uy: 0.8,
		Uz: 0.0,
	}

	q := -1.0 // electron
	mass := 1.0
	_ = units.SpeedOfLight
	bz := 2.5 // strong guide field
	dt := 0.01

	uInitialMag := math.Hypot(p.Ux, p.Uy)

	// Step for 1000 steps
	for step := 0; step < 1000; step++ {
		_, _, _, uxNew, uyNew, uzNew := pusher.PushParticle(
			&p, q, mass, dt,
			0, 0, 0, // Ex, Ey, Ez = 0
			0, 0, bz, // Bx, By, Bz
		)
		p.Ux = uxNew
		p.Uy = uyNew
		p.Uz = uzNew
	}

	uFinalMag := math.Hypot(p.Ux, p.Uy)
	relErr := math.Abs(uFinalMag-uInitialMag) / uInitialMag

	// Boris pusher guarantees machine precision conservation of magnitude
	if relErr > 1e-12 {
		t.Errorf("relativistic Boris magnetic momentum drift: relErr=%g (initial=%g, final=%g)",
			relErr, uInitialMag, uFinalMag)
	}
}

func TestRelativisticCyclotronOrbit(t *testing.T) {
	units := constants.NewNormalizedUnits(1.0, 1.0)
	pusher := particles.NewBorisPusher(units)

	c := units.SpeedOfLight
	q := 1.0
	mass := 1.0
	bz := 1.0

	// Relativistic initial velocity
	v0 := 0.6 * c
	gammaExact := 1.0 / math.Sqrt(1.0-0.6*0.6) // gamma = 1 / 0.8 = 1.25
	u0 := gammaExact * v0                      // u0 = 1.25 * 0.6 = 0.75

	// Exact relativistic cyclotron frequency and Larmor radius
	// omega_c = q * B / (gamma * m) = 1.0 * 1.0 / (1.25 * 1.0) = 0.8 rad/s
	// T_period = 2 * pi / omega_c = 2 * pi / 0.8 approx 7.85398
	// r_L = gamma * m * v0 / (q * B) = 1.25 * 1.0 * 0.6 / 1.0 = 0.75
	omegaC := (q * bz) / (gammaExact * mass)
	period := 2.0 * math.Pi / omegaC
	rLarmorExact := (gammaExact * mass * v0) / (q * bz)

	p := particles.Particle{
		X:  0.0,
		Y:  rLarmorExact, // Center of circular orbit is at (0, 0)
		Z:  0.0,
		Ux: u0,
		Uy: 0.0,
		Uz: 0.0,
	}

	nStepsPerPeriod := 200
	dt := period / float64(nStepsPerPeriod)
	numPeriods := 5

	for step := 0; step < numPeriods*nStepsPerPeriod; step++ {
		xNew, yNew, zNew, uxNew, uyNew, uzNew := pusher.PushParticle(
			&p, q, mass, dt,
			0, 0, 0,
			0, 0, bz,
		)
		p.X = xNew
		p.Y = yNew
		p.Z = zNew
		p.Ux = uxNew
		p.Uy = uyNew
		p.Uz = uzNew

		// Check radius at every step (within 2% for un-centered leapfrog finite difference)
		rCurrent := math.Hypot(p.X, p.Y)
		rErr := math.Abs(rCurrent-rLarmorExact) / rLarmorExact
		if rErr > 0.02 {
			t.Fatalf("Larmor radius deviated too much at step %d: got %g, expected %g",
				step, rCurrent, rLarmorExact)
		}
	}

	// At end of 5 complete periods, particle should be back at starting position
	finalR := math.Hypot(p.X, p.Y)
	if math.Abs(finalR-rLarmorExact)/rLarmorExact > 0.005 {
		t.Errorf("final orbital radius deviated: got %g, expected %g", finalR, rLarmorExact)
	}
}

func TestSpeciesParallelPusher(t *testing.T) {
	units := constants.NewNormalizedUnits(1.0, 1.0)
	nx := 64
	dx := 0.1
	grid := fields.NewGrid1D(nx, dx, units, fields.NewPeriodicBoundary())
	pusher := particles.NewBorisPusher(units)

	sp := particles.NewSpecies(0, "electrons", -1.0, 1.0, 100)
	for i := 0; i < 100; i++ {
		sp.AddParticle(particles.Particle{
			X:      float64(i) * 0.05,
			Ux:     0.1,
			Weight: 1.0,
		})
	}

	eKinInitial := sp.TotalKineticEnergy(units.SpeedOfLight)
	if eKinInitial <= 0 {
		t.Fatalf("expected positive kinetic energy, got %g", eKinInitial)
	}

	// Push in zero fields across 4 parallel workers
	pusher.PushSpecies1DParallel(sp, grid, 0.01, 4)

	// Since E=0 and B=0, kinetic energy must be invariant
	eKinFinal := sp.TotalKineticEnergy(units.SpeedOfLight)
	if math.Abs(eKinFinal-eKinInitial) > 1e-12 {
		t.Errorf("parallel push altered kinetic energy without fields: %g vs %g",
			eKinFinal, eKinInitial)
	}
}
