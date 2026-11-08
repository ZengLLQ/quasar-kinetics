package particles

import (
	"math"
	"sync"

	"github.com/ZengLLQ/quasar-kinetics/pkg/constants"
)

// BorisPusher integrates relativistic particle trajectories in electromagnetic fields
// using the symplectic Boris algorithm.
type BorisPusher struct {
	Units *constants.Units
}

// NewBorisPusher creates an initialized Boris particle pusher.
func NewBorisPusher(units *constants.Units) *BorisPusher {
	if units == nil {
		units = constants.NewNormalizedUnits(1.0, 1.0)
	}
	return &BorisPusher{Units: units}
}

// PushParticle advances a single particle by time step dt in local fields (E, B).
// It returns the updated (x, y, z) coordinates and (ux, uy, uz) proper velocities.
func (b *BorisPusher) PushParticle(
	p *Particle,
	q, mass, dt float64,
	ex, ey, ez float64,
	bx, by, bz float64,
) (xNew, yNew, zNew, uxNew, uyNew, uzNew float64) {
	c := b.Units.SpeedOfLight

	// Electric acceleration half-step coefficient: alpha = q * dt / (2 * m)
	alpha := (q * dt) / (2.0 * mass)

	// Step 1: Half-step electric kick: u^- = u + alpha * E
	uMinusX := p.Ux + alpha*ex
	uMinusY := p.Uy + alpha*ey
	uMinusZ := p.Uz + alpha*ez

	// Step 2: Compute gamma at intermediate state
	uMinusSq := (uMinusX*uMinusX + uMinusY*uMinusY + uMinusZ*uMinusZ) / (c * c)
	gammaMinus := math.Sqrt(1.0 + uMinusSq)

	// Step 3: Magnetic rotation vector: t = (q * dt / (2 * m * gamma)) * B
	tCoeff := alpha / gammaMinus
	tx := tCoeff * bx
	ty := tCoeff * by
	tz := tCoeff * bz

	tMagSq := tx*tx + ty*ty + tz*tz
	sCoeff := 2.0 / (1.0 + tMagSq)
	sx := sCoeff * tx
	sy := sCoeff * ty
	sz := sCoeff * tz

	// Step 4: Cross product: u' = u^- + u^- x t
	uPrimeX := uMinusX + (uMinusY*tz - uMinusZ*ty)
	uPrimeY := uMinusY + (uMinusZ*tx - uMinusX*tz)
	uPrimeZ := uMinusZ + (uMinusX*ty - uMinusY*tx)

	// Step 5: Second rotation: u^+ = u^- + u' x s
	uPlusX := uMinusX + (uPrimeY*sz - uPrimeZ*sy)
	uPlusY := uMinusY + (uPrimeZ*sx - uPrimeX*sz)
	uPlusZ := uMinusZ + (uPrimeX*sy - uPrimeY*sx)

	// Step 6: Second half-step electric kick: u^{n+1/2} = u^+ + alpha * E
	uxNew = uPlusX + alpha*ex
	uyNew = uPlusY + alpha*ey
	uzNew = uPlusZ + alpha*ez

	// Step 7: Push position using updated velocity: v = u / gamma
	uNewSq := (uxNew*uxNew + uyNew*uyNew + uzNew*uzNew) / (c * c)
	gammaNew := math.Sqrt(1.0 + uNewSq)

	vx := uxNew / gammaNew
	vy := uyNew / gammaNew
	vz := uzNew / gammaNew

	xNew = p.X + dt*vx
	yNew = p.Y + dt*vy
	zNew = p.Z + dt*vz

	return xNew, yNew, zNew, uxNew, uyNew, uzNew
}

// FieldInterpolator defines the interface required by the particle pusher to sample fields.
type FieldInterpolator interface {
	InterpolateE(x float64) (ex, ey, ez float64)
	InterpolateB(x float64) (bx, by, bz float64)
	DomainLength() float64
}

// PushSpecies1DParallel advances all particles in a species concurrently across goroutines.
func (b *BorisPusher) PushSpecies1DParallel(
	sp *Species,
	fields FieldInterpolator,
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

	lx := fields.DomainLength()
	chunkSize := (n + numWorkers - 1) / numWorkers

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

		wg.Add(1)
		go func(pStart, pEnd int) {
			defer wg.Done()
			for i := pStart; i < pEnd; i++ {
				p := &sp.Particles[i]

				// Interpolate fields at current particle position
				ex, ey, ez := fields.InterpolateE(p.X)
				bx, by, bz := fields.InterpolateB(p.X)

				// Push with Boris algorithm
				xNew, yNew, zNew, uxNew, uyNew, uzNew := b.PushParticle(
					p, sp.Charge, sp.Mass, dt,
					ex, ey, ez,
					bx, by, bz,
				)

				// Periodic boundary wrapping in x
				xWrapped := math.Mod(xNew, lx)
				if xWrapped < 0 {
					xWrapped += lx
				}

				p.X = xWrapped
				p.Y = yNew
				p.Z = zNew
				p.Ux = uxNew
				p.Uy = uyNew
				p.Uz = uzNew
			}
		}(start, end)
	}

	wg.Wait()
}
