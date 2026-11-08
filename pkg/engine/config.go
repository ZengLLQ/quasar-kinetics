package engine

import (
	"fmt"

	"github.com/ZengLLQ/quasar-kinetics/pkg/constants"
	"github.com/ZengLLQ/quasar-kinetics/pkg/fields"
)

// Config encapsulates configuration parameters for a PIC simulation run.
type Config struct {
	Nx                 int                   // Number of grid cells
	Dx                 float64               // Cell size dx
	Dt                 float64               // Time step dt (0 for auto CFL)
	CFL                float64               // Courant factor (e.g. 0.95)
	TotalSteps         int                   // Number of simulation steps
	NumWorkers         int                   // Number of worker goroutines
	Units              *constants.Units      // SI or Normalized units
	Boundary           fields.BoundaryConfig // Boundary condition settings
	DiagnosticInterval int                   // Steps between diagnostic output
}

// DefaultConfig returns recommended configuration for standard normalized simulations.
func DefaultConfig(nx int, lx float64) *Config {
	dx := lx / float64(nx)
	units := constants.NewNormalizedUnits(1.0, 1.0)
	cfl := 0.95
	dt := cfl * dx / units.SpeedOfLight

	return &Config{
		Nx:                 nx,
		Dx:                 dx,
		Dt:                 dt,
		CFL:                cfl,
		TotalSteps:         500,
		NumWorkers:         4,
		Units:              units,
		Boundary:           fields.NewPeriodicBoundary(),
		DiagnosticInterval: 10,
	}
}

// Validate checks that configuration values are physically sound and numerically stable.
func (c *Config) Validate() error {
	if c.Nx <= 2 {
		return fmt.Errorf("grid cells Nx must be > 2, got %d", c.Nx)
	}
	if c.Dx <= 0 {
		return fmt.Errorf("cell size Dx must be positive, got %g", c.Dx)
	}
	if c.CFL <= 0 || c.CFL > 1.0 {
		return fmt.Errorf("CFL parameter must be in (0, 1.0], got %g", c.CFL)
	}
	if c.Units == nil {
		c.Units = constants.NewNormalizedUnits(1.0, 1.0)
	}
	dtMax := c.CFL * c.Dx / c.Units.SpeedOfLight
	if c.Dt <= 0 {
		c.Dt = dtMax
	} else if c.Dt > c.Dx/c.Units.SpeedOfLight {
		return fmt.Errorf("time step Dt=%g violates Courant condition (c*Dt/Dx=%g > 1)",
			c.Dt, (c.Units.SpeedOfLight*c.Dt)/c.Dx)
	}
	if c.NumWorkers <= 0 {
		c.NumWorkers = 1
	}
	return nil
}
