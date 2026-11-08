package deposition

import (
	"math"
	"sync"

	"github.com/ZengLLQ/quasar-kinetics/pkg/constants"
	"github.com/ZengLLQ/quasar-kinetics/pkg/fields"
	"github.com/ZengLLQ/quasar-kinetics/pkg/particles"
)

// CurrentDepositor handles charge-conserving current and charge deposition onto Yee lattices.
type CurrentDepositor struct {
	Units *constants.Units
}

// NewCurrentDepositor creates an initialized CurrentDepositor.
func NewCurrentDepositor(units *constants.Units) *CurrentDepositor {
	if units == nil {
		units = constants.NewNormalizedUnits(1.0, 1.0)
	}
	return &CurrentDepositor{Units: units}
}

// DepositChargeCIC accumulates scalar charge density rho onto grid nodes using standard Cloud-In-Cell.
//
//	rho_j = sum_p (q_p * w_p / dx) * S(x_p - x_j)
func (d *CurrentDepositor) DepositChargeCIC(g *fields.Grid1D, sp *particles.Species) {
	dx := g.Dx
	invDx := 1.0 / dx
	coeff := (sp.Charge * invDx)
	nx := g.Nx
	lx := g.Lx

	for i := range sp.Particles {
		p := &sp.Particles[i]
		// Ensure x is wrapped within domain [0, Lx)
		xRel := math.Mod(p.X-g.X0, lx)
		if xRel < 0 {
			xRel += lx
		}

		u := xRel * invDx
		baseIdx, w0, w1 := WeightsCIC(u)

		baseIdx = ((baseIdx % nx) + nx) % nx
		nextIdx := (baseIdx + 1) % nx

		idx0 := g.Nghost + baseIdx
		idx1 := g.Nghost + nextIdx

		qEff := coeff * p.Weight
		g.Rho[idx0] += qEff * w0
		g.Rho[idx1] += qEff * w1
	}

	g.ApplyBoundaries()
}

// DepositCurrent1D accumulates currents Jx, Jy, Jz from a particle motion from x1 to x2 during time dt
// using the 1D Villasenor-Buneman / Esirkepov charge-conserving algorithm.
//
// Discrete continuity: (rho^{n+1}_j - rho^n_j)/dt + (Jx_{j+1/2} - Jx_{j-1/2})/dx = 0 is satisfied exactly.
func (d *CurrentDepositor) DepositCurrent1D(
	g *fields.Grid1D,
	x1, x2 float64,
	uxAvg, uyAvg, uzAvg float64,
	q, weight, dt float64,
) {
	dx := g.Dx
	invDx := 1.0 / dx
	invDt := 1.0 / dt
	gammaAvg := constants.RelativisticGamma(uxAvg, uyAvg, uzAvg, d.Units.SpeedOfLight)
	c := d.Units.SpeedOfLight
	_, vyAvg, vzAvg := constants.VelocityFromProperMom(uxAvg, uyAvg, uzAvg, c)

	u1 := (x1 - g.X0) * invDx
	u2 := (x2 - g.X0) * invDx

	i1, w1_0, w1_1 := WeightsCIC(u1)
	i2, w2_0, w2_1 := WeightsCIC(u2)

	minI := i1
	if i2 < minI {
		minI = i2
	}
	maxI := i1 + 1
	if i2+1 > maxI {
		maxI = i2 + 1
	}

	qwInvDt := q * weight * invDt
	_ = gammaAvg

	// Compute differences in shape function at each affected node j: D(j) = W2(j) - W1(j)
	// Then J_{j+1/2} = J_{j-1/2} - (qw / dt) * D(j)
	currJ := 0.0
	nx := g.Nx
	for j := minI; j <= maxI; j++ {
		// Evaluate W1(j) and W2(j)
		var w1j, w2j float64
		if j == i1 {
			w1j = w1_0
		} else if j == i1+1 {
			w1j = w1_1
		}

		if j == i2 {
			w2j = w2_0
		} else if j == i2+1 {
			w2j = w2_1
		}

		diff := w2j - w1j
		currJ -= qwInvDt * diff

		// Wrap index modulo physical grid
		jWrapped := ((j % nx) + nx) % nx
		idx := g.Nghost + jWrapped

		g.Jx[idx] += currJ

		// Transverse currents Jy, Jz from average velocities
		wAvg := 0.5 * (w1j + w2j)
		g.Jy[idx] += (q * weight * invDx * vyAvg) * wAvg
		g.Jz[idx] += (q * weight * invDx * vzAvg) * wAvg
	}
	g.ApplyBoundaries()
}

// LocalCurrentAccumulator holds thread-local current buffers for parallel lock-free deposition.
type LocalCurrentAccumulator struct {
	Jx  []float64
	Jy  []float64
	Jz  []float64
	Rho []float64
}

// NewLocalAccumulator allocates buffers matching grid dimensions.
func NewLocalAccumulator(totalX int) *LocalCurrentAccumulator {
	return &LocalCurrentAccumulator{
		Jx:  make([]float64, totalX),
		Jy:  make([]float64, totalX),
		Jz:  make([]float64, totalX),
		Rho: make([]float64, totalX),
	}
}

// DepositSpeciesCurrentsParallel deposits currents for an entire species using concurrent worker pools.
func (d *CurrentDepositor) DepositSpeciesCurrentsParallel(
	g *fields.Grid1D,
	sp *particles.Species,
	oldPositions []float64,
	dt float64,
	numWorkers int,
) {
	n := len(sp.Particles)
	if n == 0 {
		return
	}

	if numWorkers <= 0 {
		numWorkers = 4
	}

	chunkSize := (n + numWorkers - 1) / numWorkers
	accumulators := make([]*LocalCurrentAccumulator, numWorkers)
	for w := 0; w < numWorkers; w++ {
		accumulators[w] = NewLocalAccumulator(g.TotalX)
	}

	var wg sync.WaitGroup

	for w := 0; w < numWorkers; w++ {
		start := w * chunkSize
		end := start + chunkSize
		if start >= n {
			break
		}
		if end > n {
			end = n
		}

		acc := accumulators[w]
		wg.Add(1)
		go func(pStart, pEnd int, local *LocalCurrentAccumulator) {
			defer wg.Done()
			dx := g.Dx
			invDx := 1.0 / dx
			invDt := 1.0 / dt
			c := d.Units.SpeedOfLight

			for i := pStart; i < pEnd; i++ {
				p := &sp.Particles[i]
				x1 := oldPositions[i]
				x2 := p.X

				// Handle periodic boundary wrapping distance minimization:
				// if |x2 - x1| > Lx/2, particle wrapped
				lx := g.Lx
				dxStep := x2 - x1
				if dxStep > 0.5*lx {
					x1 += lx
				} else if dxStep < -0.5*lx {
					x1 -= lx
				}

				u1 := (x1 - g.X0) * invDx
				u2 := (x2 - g.X0) * invDx

				i1, w1_0, w1_1 := WeightsCIC(u1)
				i2, w2_0, w2_1 := WeightsCIC(u2)

				minI := i1
				if i2 < minI {
					minI = i2
				}
				maxI := i1 + 1
				if i2+1 > maxI {
					maxI = i2 + 1
				}

				qwInvDt := sp.Charge * p.Weight * invDt
				_, vyAvg, vzAvg := constants.VelocityFromProperMom(p.Ux, p.Uy, p.Uz, c)

				currJ := 0.0
				nx := g.Nx
				for j := minI; j <= maxI; j++ {
					var w1j, w2j float64
					if j == i1 {
						w1j = w1_0
					} else if j == i1+1 {
						w1j = w1_1
					}

					if j == i2 {
						w2j = w2_0
					} else if j == i2+1 {
						w2j = w2_1
					}

					diff := w2j - w1j
					currJ -= qwInvDt * diff

					jWrapped := ((j % nx) + nx) % nx
					idx := g.Nghost + jWrapped

					local.Jx[idx] += currJ

					wAvg := 0.5 * (w1j + w2j)
					local.Jy[idx] += (sp.Charge * p.Weight * invDx * vyAvg) * wAvg
					local.Jz[idx] += (sp.Charge * p.Weight * invDx * vzAvg) * wAvg
				}
			}
		}(start, end, acc)
	}

	wg.Wait()

	// Reduce all local accumulators into the global grid
	for w := 0; w < numWorkers; w++ {
		acc := accumulators[w]
		for i := 0; i < g.TotalX; i++ {
			g.Jx[i] += acc.Jx[i]
			g.Jy[i] += acc.Jy[i]
			g.Jz[i] += acc.Jz[i]
		}
	}

	g.ReduceCurrentsFromGhostCells()
}
