package fields

import (
	"math"

	"github.com/ZengLLQ/quasar-kinetics/pkg/constants"
)

// Grid1D represents a 1D3V staggered Yee electromagnetic lattice.
// Physical cells range from index [ghost] to [ghost + Nx - 1].
type Grid1D struct {
	Nx     int     // Number of physical cells
	Nghost int     // Number of ghost cells on each side (typically 2)
	TotalX int     // Total cells including ghosts: Nx + 2*Nghost
	Dx     float64 // Spatial cell size [m] or normalized
	Lx     float64 // Total physical domain length: Nx * Dx
	X0     float64 // Domain left origin

	// Electric field components (E_x is at nodes, E_y, E_z at edges)
	Ex []float64
	Ey []float64
	Ez []float64

	// Magnetic field components (B_x is constant, B_y, B_z at staggered half-nodes)
	Bx []float64
	By []float64
	Bz []float64

	// Current density components
	Jx []float64
	Jy []float64
	Jz []float64

	// Charge density and electrostatic scalar potential
	Rho []float64
	Phi []float64

	// Units and boundary configuration
	Units    *constants.Units
	Boundary BoundaryConfig
}

// NewGrid1D creates an initialized 1D3V Yee lattice.
func NewGrid1D(nx int, dx float64, units *constants.Units, bc BoundaryConfig) *Grid1D {
	ng := 2
	tot := nx + 2*ng
	if units == nil {
		units = constants.NewNormalizedUnits(1.0, 1.0)
	}

	return &Grid1D{
		Nx:       nx,
		Nghost:   ng,
		TotalX:   tot,
		Dx:       dx,
		Lx:       float64(nx) * dx,
		X0:       0.0,
		Ex:       make([]float64, tot),
		Ey:       make([]float64, tot),
		Ez:       make([]float64, tot),
		Bx:       make([]float64, tot),
		By:       make([]float64, tot),
		Bz:       make([]float64, tot),
		Jx:       make([]float64, tot),
		Jy:       make([]float64, tot),
		Jz:       make([]float64, tot),
		Rho:      make([]float64, tot),
		Phi:      make([]float64, tot),
		Units:    units,
		Boundary: bc,
	}
}

// CellCoordinate returns the physical position of cell node i.
func (g *Grid1D) CellCoordinate(i int) float64 {
	return g.X0 + float64(i-g.Nghost)*g.Dx
}

// DomainLength returns the physical length of the grid domain along x.
func (g *Grid1D) DomainLength() float64 {
	return g.Lx
}

// ResetCurrents zeroes out current density and charge arrays before deposition.
func (g *Grid1D) ResetCurrents() {
	for i := 0; i < g.TotalX; i++ {
		g.Jx[i] = 0.0
		g.Jy[i] = 0.0
		g.Jz[i] = 0.0
		g.Rho[i] = 0.0
	}
}

// InterpolateE computes electric field (Ex, Ey, Ez) at an arbitrary continuous position x
// using 1st order Cloud-in-Cell (CIC) linear weighting.
func (g *Grid1D) InterpolateE(x float64) (ex, ey, ez float64) {
	// Normalize coordinate with respect to domain
	xNorm := math.Mod(x-g.X0, g.Lx)
	if xNorm < 0 {
		xNorm += g.Lx
	}

	// Fractional index relative to physical grid
	u := xNorm / g.Dx
	i := int(math.Floor(u))
	w1 := u - float64(i)
	w0 := 1.0 - w1

	idx0 := g.Nghost + i
	idx1 := idx0 + 1
	if i >= g.Nx-1 {
		idx1 = g.Nghost // Wrap for periodic
	}

	ex = w0*g.Ex[idx0] + w1*g.Ex[idx1]
	ey = w0*g.Ey[idx0] + w1*g.Ey[idx1]
	ez = w0*g.Ez[idx0] + w1*g.Ez[idx1]
	return ex, ey, ez
}

// InterpolateB computes magnetic field (Bx, By, Bz) at continuous position x
// taking into account spatial staggering of half-nodes.
func (g *Grid1D) InterpolateB(x float64) (bx, by, bz float64) {
	xNorm := math.Mod(x-g.X0, g.Lx)
	if xNorm < 0 {
		xNorm += g.Lx
	}

	// Magnetic fields By, Bz are defined on staggered half-nodes x_{i+1/2} = (i + 0.5) * dx
	uStaggered := (xNorm / g.Dx) - 0.5
	i := int(math.Floor(uStaggered))
	w1 := uStaggered - float64(i)
	w0 := 1.0 - w1

	idx0 := g.Nghost + i
	idx1 := idx0 + 1
	if i < 0 {
		idx0 = g.Nghost + g.Nx - 1
		idx1 = g.Nghost
	} else if i >= g.Nx-1 {
		idx1 = g.Nghost
	}

	bx = g.Bx[g.Nghost] // Bx is constant 1D background guide field
	by = w0*g.By[idx0] + w1*g.By[idx1]
	bz = w0*g.Bz[idx0] + w1*g.Bz[idx1]
	return bx, by, bz
}

// ElectricEnergy calculates total electric field energy in the physical domain:
//
//	U_E = (eps_0 / 2) * sum( |E|^2 ) * dx
func (g *Grid1D) ElectricEnergy() float64 {
	sum := 0.0
	start := g.Nghost
	end := g.Nghost + g.Nx
	for i := start; i < end; i++ {
		sum += g.Ex[i]*g.Ex[i] + g.Ey[i]*g.Ey[i] + g.Ez[i]*g.Ez[i]
	}
	return 0.5 * g.Units.Permittivity * sum * g.Dx
}

// MagneticEnergy calculates total magnetic field energy in the physical domain:
//
//	U_B = (1 / (2 * mu_0)) * sum( |B|^2 ) * dx
func (g *Grid1D) MagneticEnergy() float64 {
	sum := 0.0
	start := g.Nghost
	end := g.Nghost + g.Nx
	for i := start; i < end; i++ {
		sum += g.Bx[i]*g.Bx[i] + g.By[i]*g.By[i] + g.Bz[i]*g.Bz[i]
	}
	return (0.5 / g.Units.Permeability) * sum * g.Dx
}

// ApplyBoundaries synchronizes ghost cells according to the boundary configuration.
func (g *Grid1D) ApplyBoundaries() {
	ng := g.Nghost
	nx := g.Nx

	if g.Boundary.XLeft == BoundaryPeriodic && g.Boundary.XRight == BoundaryPeriodic {
		// Periodic boundary wrap
		for k := 0; k < ng; k++ {
			// Left ghost cells copied from right interior
			g.Ex[k] = g.Ex[nx+k]
			g.Ey[k] = g.Ey[nx+k]
			g.Ez[k] = g.Ez[nx+k]
			g.Bx[k] = g.Bx[nx+k]
			g.By[k] = g.By[nx+k]
			g.Bz[k] = g.Bz[nx+k]
			g.Jx[k] = g.Jx[nx+k]
			g.Jy[k] = g.Jy[nx+k]
			g.Jz[k] = g.Jz[nx+k]
			g.Rho[k] = g.Rho[nx+k]

			// Right ghost cells copied from left interior
			g.Ex[nx+ng+k] = g.Ex[ng+k]
			g.Ey[nx+ng+k] = g.Ey[ng+k]
			g.Ez[nx+ng+k] = g.Ez[ng+k]
			g.Bx[nx+ng+k] = g.Bx[ng+k]
			g.By[nx+ng+k] = g.By[ng+k]
			g.Bz[nx+ng+k] = g.Bz[ng+k]
			g.Jx[nx+ng+k] = g.Jx[ng+k]
			g.Jy[nx+ng+k] = g.Jy[ng+k]
			g.Jz[nx+ng+k] = g.Jz[ng+k]
			g.Rho[nx+ng+k] = g.Rho[ng+k]
		}
	} else if g.Boundary.XLeft == BoundaryConducting {
		// Conducting boundary: tangential E vanishes at boundary
		g.Ey[ng] = 0.0
		g.Ez[ng] = 0.0
		g.Ey[ng+nx-1] = 0.0
		g.Ez[ng+nx-1] = 0.0
	}
}

// ReduceCurrentsFromGhostCells folds ghost cell deposited currents back into physical cells for periodic grids.
func (g *Grid1D) ReduceCurrentsFromGhostCells() {
	if g.Boundary.XLeft != BoundaryPeriodic {
		return
	}
	ng := g.Nghost
	nx := g.Nx

	for k := 0; k < ng; k++ {
		// Ghost cells beyond right boundary fold into left physical cells
		g.Jx[ng+k] += g.Jx[nx+ng+k]
		g.Jy[ng+k] += g.Jy[nx+ng+k]
		g.Jz[ng+k] += g.Jz[nx+ng+k]
		g.Rho[ng+k] += g.Rho[nx+ng+k]

		// Ghost cells before left boundary fold into right physical cells
		g.Jx[nx+k] += g.Jx[k]
		g.Jy[nx+k] += g.Jy[k]
		g.Jz[nx+k] += g.Jz[k]
		g.Rho[nx+k] += g.Rho[k]
	}

	// Update boundaries again to ensure periodicity
	g.ApplyBoundaries()
}
