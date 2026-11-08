package deposition_test

import (
	"math"
	"testing"

	"github.com/ZengLLQ/quasar-kinetics/pkg/constants"
	"github.com/ZengLLQ/quasar-kinetics/pkg/deposition"
	"github.com/ZengLLQ/quasar-kinetics/pkg/fields"
	"github.com/ZengLLQ/quasar-kinetics/pkg/particles"
)

func TestEsirkepovExactChargeConservation(t *testing.T) {
	units := constants.NewNormalizedUnits(1.0, 1.0)
	nx := 32
	dx := 0.2
	grid := fields.NewGrid1D(nx, dx, units, fields.NewPeriodicBoundary())
	depositor := deposition.NewCurrentDepositor(units)

	q := -1.0
	weight := 10.0
	dt := 0.05

	// Test multiple particle displacements within and across cell boundaries
	testCases := []struct {
		x1, x2 float64
	}{
		{x1: 1.05, x2: 1.15}, // Small motion within same cell
		{x1: 1.10, x2: 1.35}, // Motion crossing cell boundary
		{x1: 2.35, x2: 2.10}, // Backward motion
		{x1: 0.50, x2: 0.90}, // Spanning two cell widths
	}

	for caseIdx, tc := range testCases {
		grid.ResetCurrents()

		// 1. Deposit charge at initial position x1
		sp1 := particles.NewSpecies(0, "e1", q, 1.0, 1)
		sp1.AddParticle(particles.Particle{X: tc.x1, Weight: weight})
		depositor.DepositChargeCIC(grid, sp1)
		rho1 := make([]float64, grid.TotalX)
		copy(rho1, grid.Rho)

		// 2. Deposit charge at final position x2
		grid.ResetCurrents()
		sp2 := particles.NewSpecies(0, "e2", q, 1.0, 1)
		sp2.AddParticle(particles.Particle{X: tc.x2, Weight: weight})
		depositor.DepositChargeCIC(grid, sp2)
		rho2 := make([]float64, grid.TotalX)
		copy(rho2, grid.Rho)

		// 3. Deposit current from trajectory x1 -> x2
		grid.ResetCurrents()
		depositor.DepositCurrent1D(grid, tc.x1, tc.x2, 0.5, 0, 0, q, weight, dt)

		// 4. Verify discrete continuity equation at every node j:
		// d_rho/dt + div(J) = 0
		// (rho2[j] - rho1[j])/dt + (Jx[j] - Jx[j-1])/dx = 0
		var maxContinuityResidual float64
		for j := grid.Nghost; j < grid.Nghost+grid.Nx; j++ {
			dRhoDt := (rho2[j] - rho1[j]) / dt
			divJ := (grid.Jx[j] - grid.Jx[j-1]) / dx

			residual := math.Abs(dRhoDt + divJ)
			if residual > maxContinuityResidual {
				maxContinuityResidual = residual
			}
		}

		if maxContinuityResidual > 1e-12 {
			t.Errorf("case %d (x1=%g, x2=%g): discrete continuity violated! residual=%g",
				caseIdx, tc.x1, tc.x2, maxContinuityResidual)
		}
	}
}

func TestDepositTotalChargeConservation(t *testing.T) {
	units := constants.NewNormalizedUnits(1.0, 1.0)
	nx := 64
	dx := 0.1
	grid := fields.NewGrid1D(nx, dx, units, fields.NewPeriodicBoundary())
	depositor := deposition.NewCurrentDepositor(units)

	sp := particles.NewSpecies(0, "ions", 1.0, 1836.0, 50)
	var expectedTotalCharge float64
	for i := 0; i < 50; i++ {
		w := float64(i + 1)
		sp.AddParticle(particles.Particle{
			X:      float64(i) * 0.12,
			Weight: w,
		})
		expectedTotalCharge += w * sp.Charge
	}

	depositor.DepositChargeCIC(grid, sp)

	// Integrated charge: sum(rho[i]) * dx
	var integratedCharge float64
	for i := grid.Nghost; i < grid.Nghost+grid.Nx; i++ {
		integratedCharge += grid.Rho[i] * grid.Dx
	}

	relErr := math.Abs(integratedCharge-expectedTotalCharge) / expectedTotalCharge
	if relErr > 1e-12 {
		t.Errorf("total charge mismatch: integrated=%g, expected=%g, relErr=%g",
			integratedCharge, expectedTotalCharge, relErr)
	}
}
