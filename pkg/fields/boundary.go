package fields

// BoundaryType defines the boundary conditions applied at domain edges.
type BoundaryType int

const (
	// BoundaryPeriodic wraps fields periodically across domain boundaries.
	BoundaryPeriodic BoundaryType = iota
	// BoundaryConducting sets tangential E and normal B to zero (Perfect Electric Conductor).
	BoundaryConducting
	// BoundaryAbsorbing applies Silver-Muller radiation absorbing condition.
	BoundaryAbsorbing
	// BoundaryOpen leaves ghost cells extrapolated without reflection suppression.
	BoundaryOpen
)

// BoundaryConfig specifies boundary behavior along each coordinate axis.
type BoundaryConfig struct {
	XLeft  BoundaryType
	XRight BoundaryType
	YLeft  BoundaryType
	YRight BoundaryType
}

// NewPeriodicBoundary returns a configuration with periodic boundaries along all axes.
func NewPeriodicBoundary() BoundaryConfig {
	return BoundaryConfig{
		XLeft:  BoundaryPeriodic,
		XRight: BoundaryPeriodic,
		YLeft:  BoundaryPeriodic,
		YRight: BoundaryPeriodic,
	}
}

// NewConductingBoundary returns PEC boundary configuration.
func NewConductingBoundary() BoundaryConfig {
	return BoundaryConfig{
		XLeft:  BoundaryConducting,
		XRight: BoundaryConducting,
		YLeft:  BoundaryConducting,
		YRight: BoundaryConducting,
	}
}
