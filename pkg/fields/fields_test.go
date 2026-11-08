package fields_test

import (
	"math"
	"testing"

	"github.com/ZengLLQ/quasar-kinetics/pkg/constants"
	"github.com/ZengLLQ/quasar-kinetics/pkg/fields"
)

func TestGrid1DInitialization(t *testing.T) {
	units := constants.NewNormalizedUnits(1.0, 1.0)
	nx := 64
	dx := 0.1
	bc := fields.NewPeriodicBoundary()

	grid := fields.NewGrid1D(nx, dx, units, bc)

	if grid.Nx != nx {
		t.Errorf("expected Nx=%d, got %d", nx, grid.Nx)
	}
	if grid.Lx != float64(nx)*dx {
		t.Errorf("expected Lx=%g, got %g", float64(nx)*dx, grid.Lx)
	}

	// Verify coordinates
	coord0 := grid.CellCoordinate(grid.Nghost)
	if coord0 != 0.0 {
		t.Errorf("expected origin coordinate 0.0, got %g", coord0)
	}
}

func TestInterpolationCIC(t *testing.T) {
	units := constants.NewNormalizedUnits(1.0, 1.0)
	nx := 32
	dx := 1.0
	grid := fields.NewGrid1D(nx, dx, units, fields.NewPeriodicBoundary())

	// Set linear ramp in Ex: Ex[i] = i * 2.0
	for i := 0; i < grid.TotalX; i++ {
		grid.Ex[i] = float64(i) * 2.0
	}
	grid.ApplyBoundaries()

	// Sample midpoint between node 3 and 4 (relative to interior)
	// Physical x = 3.5
	ex, _, _ := grid.InterpolateE(3.5)
	expected := float64(grid.Nghost+3)*2.0*0.5 + float64(grid.Nghost+4)*2.0*0.5
	if math.Abs(ex-expected) > 1e-12 {
		t.Errorf("expected interpolated Ex=%g, got %g", expected, ex)
	}
}

func TestPoissonSolverAnalytical(t *testing.T) {
	units := constants.NewNormalizedUnits(1.0, 1.0)
	nx := 128
	lx := 2.0 * math.Pi
	dx := lx / float64(nx)

	grid := fields.NewGrid1D(nx, dx, units, fields.NewPeriodicBoundary())
	solver := fields.NewFDTDSolver(0.95)

	// Single harmonic charge density: rho(x) = rho0 * sin(k * x)
	// d^2 phi / dx^2 = - rho / eps0 => phi(x) = (rho0 / (k^2 * eps0)) * sin(k * x)
	// E_x(x) = - d phi / dx = - (rho0 / (k * eps0)) * cos(k * x)
	rho0 := 1.0
	kMode := 1.0 // wavenumber k = 2 * pi / L = 1.0
	ng := grid.Nghost

	for i := 0; i < nx; i++ {
		x := float64(i) * dx
		grid.Rho[ng+i] = rho0 * math.Sin(kMode*x)
	}
	grid.ApplyBoundaries()

	solver.SolvePoisson(grid)

	// Check accuracy of potential and electric field against analytical solution
	var maxErrPhi, maxErrEx float64
	for i := 0; i < nx; i++ {
		x := float64(i) * dx
		phiExact := (rho0 / (kMode * kMode * units.Permittivity)) * math.Sin(kMode*x)
		exExact := -(rho0 / (kMode * units.Permittivity)) * math.Cos(kMode*x)

		errPhi := math.Abs(grid.Phi[ng+i] - phiExact)
		errEx := math.Abs(grid.Ex[ng+i] - exExact)

		if errPhi > maxErrPhi {
			maxErrPhi = errPhi
		}
		if errEx > maxErrEx {
			maxErrEx = errEx
		}
	}

	// For 128 points, discrete Laplacian eigenvalue exact inversion achieves < 2e-3 accuracy for central difference
	if maxErrPhi > 1e-2 {
		t.Errorf("Poisson potential error too large: %g", maxErrPhi)
	}
	if maxErrEx > 0.05 {
		t.Errorf("Poisson electric field error too large: %g", maxErrEx)
	}
}

func TestVacuumWavePropagationEnergyConservation(t *testing.T) {
	units := constants.NewNormalizedUnits(1.0, 1.0)
	nx := 100
	dx := 0.1
	lx := float64(nx) * dx
	grid := fields.NewGrid1D(nx, dx, units, fields.NewPeriodicBoundary())

	solver := fields.NewFDTDSolver(0.5) // CFL = 0.5 for high accuracy
	dt := solver.MaxStableTimeStep(grid)

	// Initialize right-propagating EM wave:
	// Ez(x) = E0 * sin(k * x)
	// By(x) = (E0 / c) * sin(k * x)
	k := 2.0 * math.Pi / lx
	e0 := 1.0
	ng := grid.Nghost
	for i := 0; i < nx; i++ {
		x := float64(i) * dx
		grid.Ez[ng+i] = e0 * math.Sin(k*x)
		// Staggered half-node for By: x_{i+1/2} = (i + 0.5) * dx
		grid.By[ng+i] = (e0 / units.SpeedOfLight) * math.Sin(k*(x+0.5*dx))
	}
	grid.ApplyBoundaries()

	initialEnergy := grid.ElectricEnergy() + grid.MagneticEnergy()
	if initialEnergy <= 0 {
		t.Fatalf("expected positive initial EM energy, got %g", initialEnergy)
	}

	// Propagate wave for one full period T = lx / c
	periodSteps := int(math.Round(lx / (units.SpeedOfLight * dt)))
	for step := 0; step < periodSteps; step++ {
		// Leapfrog FDTD: push B by dt/2, push E by dt, push B by dt/2
		solver.PushB(grid, 0.5*dt)
		solver.PushE(grid, dt)
		solver.PushB(grid, 0.5*dt)
	}

	finalEnergy := grid.ElectricEnergy() + grid.MagneticEnergy()
	relEnergyDrift := math.Abs(finalEnergy-initialEnergy) / initialEnergy

	if relEnergyDrift > 1e-4 {
		t.Errorf("vacuum EM energy conservation drift too high: %g (initial=%g, final=%g)",
			relEnergyDrift, initialEnergy, finalEnergy)
	}
}
