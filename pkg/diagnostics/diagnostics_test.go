package diagnostics_test

import (
	"math"
	"strings"
	"testing"

	"github.com/ZengLLQ/quasar-kinetics/pkg/constants"
	"github.com/ZengLLQ/quasar-kinetics/pkg/diagnostics"
	"github.com/ZengLLQ/quasar-kinetics/pkg/engine"
	"github.com/ZengLLQ/quasar-kinetics/pkg/fields"
	"github.com/ZengLLQ/quasar-kinetics/pkg/particles"
)

func TestPhaseSpaceBinningAndASCII(t *testing.T) {
	ps := diagnostics.NewPhaseSpace2D(20, 20, 0, 10.0, -1.0, 1.0)

	sp := particles.NewSpecies(0, "e", -1.0, 1.0, 10)
	sp.AddParticle(particles.Particle{X: 5.0, Ux: 0.5, Weight: 1.0})
	sp.AddParticle(particles.Particle{X: 2.0, Ux: -0.5, Weight: 2.0})

	ps.BinSpecies(sp, 1.0)

	// ASCII output
	ascii := ps.ASCIIHeatmap(20, 10)
	if !strings.Contains(ascii, "Phase Space f(x, vx)") {
		t.Errorf("expected header in ASCII output, got: %s", ascii)
	}
}

func TestSpectralAnalyzerSinusoid(t *testing.T) {
	nx := 64
	lx := 2.0 * math.Pi
	dx := lx / float64(nx)
	units := constants.NewNormalizedUnits(1.0, 1.0)
	grid := fields.NewGrid1D(nx, dx, units, fields.NewPeriodicBoundary())

	// Set mode m=3: Ex(x) = 2.5 * cos(3 * x)
	mTarget := 3
	aTarget := 2.5
	kTarget := float64(mTarget) * (2.0 * math.Pi / lx)

	for i := 0; i < nx; i++ {
		x := float64(i) * dx
		grid.Ex[grid.Nghost+i] = aTarget * math.Cos(kTarget*x)
	}

	sa := diagnostics.NewSpectralAnalyzer()
	modes := sa.AnalyzeExModes(grid)

	if len(modes) <= mTarget {
		t.Fatalf("not enough modes analyzed")
	}

	dom := sa.DominantWaveMode(grid)
	if dom.ModeNumber != mTarget {
		t.Errorf("expected dominant mode %d, got %d", mTarget, dom.ModeNumber)
	}

	if math.Abs(dom.Amplitude-aTarget) > 1e-10 {
		t.Errorf("expected amplitude %g, got %g", aTarget, dom.Amplitude)
	}
}

func TestEnergyReportCSVAndSummary(t *testing.T) {
	snaps := []engine.EnergySnapshot{
		{Step: 0, Time: 0.0, KineticEnergy: 10.0, ElectricEnergy: 0.0, MagneticEnergy: 0.0, TotalEnergy: 10.0, RelDrift: 0.0},
		{Step: 10, Time: 0.1, KineticEnergy: 9.8, ElectricEnergy: 0.2, MagneticEnergy: 0.0, TotalEnergy: 10.0, RelDrift: 0.001},
	}

	rep := diagnostics.NewEnergyReport(snaps)
	csv := rep.ToCSV()
	if !strings.HasPrefix(csv, "step,time,kinetic_energy") {
		t.Errorf("unexpected CSV header: %s", csv)
	}

	sum := rep.SummaryString()
	if !strings.Contains(sum, "E_kin") {
		t.Errorf("unexpected summary string: %s", sum)
	}
}
