package particles

import (
	"math"

	"github.com/ZengLLQ/quasar-kinetics/pkg/constants"
)

// Particle represents an individual computational macroparticle in phase space.
// Momenta are stored as proper velocities u = gamma * v.
type Particle struct {
	X  float64 // Spatial position x [m or normalized]
	Y  float64 // Spatial position y
	Z  float64 // Spatial position z
	Ux float64 // Proper velocity u_x = gamma * v_x
	Uy float64 // Proper velocity u_y = gamma * v_y
	Uz float64 // Proper velocity u_z = gamma * v_z

	Weight    float64 // Number of physical particles represented
	SpeciesID int     // Species identifier index
	Tag       uint64  // Unique particle tracking tag
}

// Gamma computes the Lorentz factor for this particle.
func (p *Particle) Gamma(c float64) float64 {
	return constants.RelativisticGamma(p.Ux, p.Uy, p.Uz, c)
}

// Velocity returns the 3-velocity (vx, vy, vz).
func (p *Particle) Velocity(c float64) (vx, vy, vz float64) {
	return constants.VelocityFromProperMom(p.Ux, p.Uy, p.Uz, c)
}

// KineticEnergy computes relativistic kinetic energy: E_k = m * c^2 * (gamma - 1).
func (p *Particle) KineticEnergy(mass, c float64) float64 {
	gamma := p.Gamma(c)
	return p.Weight * mass * c * c * (gamma - 1.0)
}

// Species represents a collection of macroparticles with identical charge and mass.
type Species struct {
	ID        int
	Name      string
	Charge    float64 // Particle charge q [C or normalized]
	Mass      float64 // Particle rest mass m [kg or normalized]
	Particles []Particle
}

// NewSpecies creates a new particle species with pre-allocated capacity.
func NewSpecies(id int, name string, charge, mass float64, capacity int) *Species {
	return &Species{
		ID:        id,
		Name:      name,
		Charge:    charge,
		Mass:      mass,
		Particles: make([]Particle, 0, capacity),
	}
}

// AddParticle appends a particle to the species.
func (s *Species) AddParticle(p Particle) {
	p.SpeciesID = s.ID
	s.Particles = append(s.Particles, p)
}

// TotalKineticEnergy computes the sum of kinetic energies across all macroparticles.
func (s *Species) TotalKineticEnergy(c float64) float64 {
	var total float64
	for i := range s.Particles {
		total += s.Particles[i].KineticEnergy(s.Mass, c)
	}
	return total
}

// TotalCharge computes the total physical charge represented by this species.
func (s *Species) TotalCharge() float64 {
	var total float64
	for i := range s.Particles {
		total += s.Particles[i].Weight * s.Charge
	}
	return total
}

// VelocityMoments computes mean velocity and thermal spread (temperature proxy).
func (s *Species) VelocityMoments(c float64) (meanVx, meanVy, meanVz, vThermal float64) {
	n := len(s.Particles)
	if n == 0 {
		return 0, 0, 0, 0
	}

	var sumVx, sumVy, sumVz float64
	for i := range s.Particles {
		vx, vy, vz := s.Particles[i].Velocity(c)
		sumVx += vx
		sumVy += vy
		sumVz += vz
	}

	meanVx = sumVx / float64(n)
	meanVy = sumVy / float64(n)
	meanVz = sumVz / float64(n)

	var varSum float64
	for i := range s.Particles {
		vx, vy, vz := s.Particles[i].Velocity(c)
		dvx := vx - meanVx
		dvy := vy - meanVy
		dvz := vz - meanVz
		varSum += dvx*dvx + dvy*dvy + dvz*dvz
	}

	vThermal = math.Sqrt(varSum / (3.0 * float64(n)))
	return meanVx, meanVy, meanVz, vThermal
}
