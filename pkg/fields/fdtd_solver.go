package fields

import (
	"fmt"
	"math"
)

// FDTDSolver implements the Finite-Difference Time-Domain Maxwell and Poisson equations on Yee grids.
type FDTDSolver struct {
	CFL float64 // Courant-Friedrichs-Lewy safety parameter (e.g. 0.95)
}

// NewFDTDSolver creates an FDTD solver instance.
func NewFDTDSolver(cfl float64) *FDTDSolver {
	if cfl <= 0 || cfl > 1.0 {
		cfl = 0.95
	}
	return &FDTDSolver{CFL: cfl}
}

// MaxStableTimeStep computes the maximum stable time step dt according to the CFL criterion:
//
//	dt_max = CFL * dx / c
func (s *FDTDSolver) MaxStableTimeStep(g *Grid1D) float64 {
	return s.CFL * g.Dx / g.Units.SpeedOfLight
}

// CheckCFL validates that dt satisfies the CFL condition.
func (s *FDTDSolver) CheckCFL(g *Grid1D, dt float64) error {
	dtMax := g.Dx / g.Units.SpeedOfLight
	if dt > dtMax {
		return fmt.Errorf("FDTD CFL condition violated: dt=%g exceeds dt_max=%g (dx=%g, c=%g)",
			dt, dtMax, g.Dx, g.Units.SpeedOfLight)
	}
	return nil
}

// PushB advances the magnetic field by time interval dt using Faraday's law:
//
//	curl(E) = -dB/dt
//	dB_y/dt =  dE_z/dx
//	dB_z/dt = -dE_y/dx
func (s *FDTDSolver) PushB(g *Grid1D, dt float64) {
	cfl := dt / g.Dx
	ng := g.Nghost
	nx := g.Nx

	// Interior physical half-nodes
	for i := ng; i < ng+nx; i++ {
		// B_y(i+1/2) updated by difference in E_z
		dEz := g.Ez[i+1] - g.Ez[i]
		g.By[i] += cfl * dEz

		// B_z(i+1/2) updated by difference in E_y
		dEy := g.Ey[i+1] - g.Ey[i]
		g.Bz[i] -= cfl * dEy
	}

	g.ApplyBoundaries()
}

// PushE advances the electric field by time interval dt using Ampere-Maxwell law:
//
//	dE/dt = c^2 curl(B) - J / eps_0
//	dE_x/dt = -J_x / eps_0
//	dE_y/dt = -c^2 dB_z/dx - J_y / eps_0
//	dE_z/dt =  c^2 dB_y/dx - J_z / eps_0
func (s *FDTDSolver) PushE(g *Grid1D, dt float64) {
	c2 := g.Units.SpeedOfLight * g.Units.SpeedOfLight
	cflCoeff := c2 * dt / g.Dx
	invEps := dt / g.Units.Permittivity
	ng := g.Nghost
	nx := g.Nx

	for i := ng; i < ng+nx; i++ {
		// Electrostatic component updated by longitudinal current
		g.Ex[i] -= invEps * g.Jx[i]

		// Transverse components updated by magnetic curl and currents
		dBz := g.Bz[i] - g.Bz[i-1]
		dBy := g.By[i] - g.By[i-1]

		g.Ey[i] -= cflCoeff*dBz + invEps*g.Jy[i]
		g.Ez[i] += cflCoeff*dBy - invEps*g.Jz[i]
	}

	g.ApplyBoundaries()
}

// SolvePoisson solves the 1D Poisson equation for electrostatic potential phi:
//
//	d^2 phi / dx^2 = - (rho - rho_mean) / epsilon_0
//
// using Discrete Fourier Transform (DFT/FFT) spectral decomposition with periodic boundaries,
// and derives the electrostatic field:
//
//	E_x = - d phi / dx
func (s *FDTDSolver) SolvePoisson(g *Grid1D) {
	nx := g.Nx
	ng := g.Nghost

	// Compute average net charge for neutral plasma background subtraction
	var totalCharge float64
	for i := 0; i < nx; i++ {
		totalCharge += g.Rho[ng+i]
	}
	rhoMean := totalCharge / float64(nx)

	// Subtract background charge to ensure solvability of periodic Poisson equation
	rhs := make([]float64, nx)
	for i := 0; i < nx; i++ {
		rhs[i] = -(g.Rho[ng+i] - rhoMean) / g.Units.Permittivity
	}

	// Solve via Discrete Fourier Transform
	// d^2 phi / dx^2 on discrete grid has eigenvalue:
	// lambda_k = 2 * (cos(2 * pi * k / N) - 1) / dx^2
	dx := g.Dx
	phiInterior := make([]float64, nx)

	// Forward DFT
	cosTable := make([][]float64, nx)
	sinTable := make([][]float64, nx)
	for k := 0; k < nx; k++ {
		cosTable[k] = make([]float64, nx)
		sinTable[k] = make([]float64, nx)
		for j := 0; j < nx; j++ {
			theta := 2.0 * math.Pi * float64(k*j) / float64(nx)
			cosTable[k][j] = math.Cos(theta)
			sinTable[k][j] = math.Sin(theta)
		}
	}

	realDFT := make([]float64, nx)
	imagDFT := make([]float64, nx)

	for k := 0; k < nx; k++ {
		var rSum, iSum float64
		for j := 0; j < nx; j++ {
			rSum += rhs[j] * cosTable[k][j]
			iSum -= rhs[j] * sinTable[k][j]
		}
		realDFT[k] = rSum
		imagDFT[k] = iSum
	}

	// Invert in Fourier space
	realPhiHat := make([]float64, nx)
	imagPhiHat := make([]float64, nx)
	for k := 1; k < nx; k++ {
		theta := 2.0 * math.Pi * float64(k) / float64(nx)
		denom := 2.0 * (math.Cos(theta) - 1.0) / (dx * dx)
		if math.Abs(denom) > 1e-15 {
			realPhiHat[k] = realDFT[k] / denom
			imagPhiHat[k] = imagDFT[k] / denom
		}
	}
	// k=0 mode is set to zero (gauge condition)
	realPhiHat[0] = 0.0
	imagPhiHat[0] = 0.0

	// Inverse DFT
	invN := 1.0 / float64(nx)
	for j := 0; j < nx; j++ {
		var sum float64
		for k := 0; k < nx; k++ {
			sum += realPhiHat[k]*cosTable[k][j] - imagPhiHat[k]*sinTable[k][j]
		}
		phiInterior[j] = sum * invN
	}

	// Copy into grid and compute E_x = - d phi / dx
	for j := 0; j < nx; j++ {
		g.Phi[ng+j] = phiInterior[j]
	}

	// Periodic boundary for phi
	for k := 0; k < ng; k++ {
		g.Phi[k] = g.Phi[nx+k]
		g.Phi[nx+ng+k] = g.Phi[ng+k]
	}

	// Central difference for Ex: E_x(i) = - (phi(i+1) - phi(i-1)) / (2 * dx)
	for i := ng; i < ng+nx; i++ {
		g.Ex[i] = -(g.Phi[i+1] - g.Phi[i-1]) / (2.0 * dx)
	}

	g.ApplyBoundaries()
}
