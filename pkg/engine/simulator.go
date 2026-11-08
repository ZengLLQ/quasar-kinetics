package engine

import (
	"math"

	"github.com/ZengLLQ/quasar-kinetics/pkg/deposition"
	"github.com/ZengLLQ/quasar-kinetics/pkg/fields"
	"github.com/ZengLLQ/quasar-kinetics/pkg/particles"
)

// EnergySnapshot records instantaneous energy metrics of the simulation.
type EnergySnapshot struct {
	Step           int     `json:"step"`
	Time           float64 `json:"time"`
	KineticEnergy  float64 `json:"kinetic_energy"`
	ElectricEnergy float64 `json:"electric_energy"`
	MagneticEnergy float64 `json:"magnetic_energy"`
	TotalEnergy    float64 `json:"total_energy"`
	RelDrift       float64 `json:"rel_drift"`
}

// Simulator coordinates the Yee grid, Maxwell-FDTD field solvers, Boris pushers,
// and Esirkepov current deposition into a self-consistent kinetic time loop.
type Simulator struct {
	Config    *Config
	Grid      *fields.Grid1D
	Solver    *fields.FDTDSolver
	Pusher    *particles.BorisPusher
	Depositor *deposition.CurrentDepositor
	Species   []*particles.Species

	StepCount     int
	Time          float64
	InitialEnergy float64
	History       []EnergySnapshot

	// Temporary buffers for particle positions before push
	oldPosBuffer [][]float64
}

// NewSimulator constructs and initializes a new Particle-in-Cell simulator.
func NewSimulator(cfg *Config) (*Simulator, error) {
	if err := cfg.Validate(); err != nil {
		return nil, err
	}

	grid := fields.NewGrid1D(cfg.Nx, cfg.Dx, cfg.Units, cfg.Boundary)
	solver := fields.NewFDTDSolver(cfg.CFL)
	pusher := particles.NewBorisPusher(cfg.Units)
	depositor := deposition.NewCurrentDepositor(cfg.Units)

	return &Simulator{
		Config:       cfg,
		Grid:         grid,
		Solver:       solver,
		Pusher:       pusher,
		Depositor:    depositor,
		Species:      make([]*particles.Species, 0, 4),
		StepCount:    0,
		Time:         0.0,
		History:      make([]EnergySnapshot, 0, cfg.TotalSteps/cfg.DiagnosticInterval+1),
		oldPosBuffer: make([][]float64, 0, 4),
	}, nil
}

// AddSpecies registers a particle species into the simulation engine.
func (s *Simulator) AddSpecies(sp *particles.Species) {
	s.Species = append(s.Species, sp)
	s.oldPosBuffer = append(s.oldPosBuffer, make([]float64, 0, len(sp.Particles)))
}

// InitializeElectrostatic deposits initial particle charges and solves Poisson's equation
// to establish initial self-consistent electrostatic fields E_x(x).
func (s *Simulator) InitializeElectrostatic() {
	s.Grid.ResetCurrents()
	for _, sp := range s.Species {
		s.Depositor.DepositChargeCIC(s.Grid, sp)
	}
	s.Solver.SolvePoisson(s.Grid)
	s.RecordEnergy()
}

// TotalEnergy evaluates the combined particle kinetic and electromagnetic field energy.
func (s *Simulator) TotalEnergy() (eKin, eElec, eMag, eTot float64) {
	c := s.Config.Units.SpeedOfLight
	for _, sp := range s.Species {
		eKin += sp.TotalKineticEnergy(c)
	}
	eElec = s.Grid.ElectricEnergy()
	eMag = s.Grid.MagneticEnergy()
	eTot = eKin + eElec + eMag
	return eKin, eElec, eMag, eTot
}

// RecordEnergy logs an instantaneous energy snapshot into history.
func (s *Simulator) RecordEnergy() EnergySnapshot {
	eKin, eElec, eMag, eTot := s.TotalEnergy()

	if s.InitialEnergy == 0 && eTot > 0 {
		s.InitialEnergy = eTot
	}

	var drift float64
	if s.InitialEnergy > 0 {
		drift = math.Abs(eTot-s.InitialEnergy) / s.InitialEnergy
	}

	snap := EnergySnapshot{
		Step:           s.StepCount,
		Time:           s.Time,
		KineticEnergy:  eKin,
		ElectricEnergy: eElec,
		MagneticEnergy: eMag,
		TotalEnergy:    eTot,
		RelDrift:       drift,
	}
	s.History = append(s.History, snap)
	return snap
}

// Step advances the full simulation by one time interval dt using the leapfrog PIC cycle.
func (s *Simulator) Step() {
	dt := s.Config.Dt
	halfDt := 0.5 * dt
	workers := s.Config.NumWorkers

	// Step 1: Magnetic field half-push: B^{n-1/2} -> B^n
	s.Solver.PushB(s.Grid, halfDt)

	// Step 2: Store previous positions and advance particles with Boris pusher
	for spIdx, sp := range s.Species {
		n := len(sp.Particles)
		// Ensure old positions buffer capacity
		if len(s.oldPosBuffer[spIdx]) < n {
			s.oldPosBuffer[spIdx] = make([]float64, n)
		} else {
			s.oldPosBuffer[spIdx] = s.oldPosBuffer[spIdx][:n]
		}

		for i := 0; i < n; i++ {
			s.oldPosBuffer[spIdx][i] = sp.Particles[i].X
		}

		// Advance particle velocities and positions
		s.Pusher.PushSpecies1DParallel(sp, s.Grid, dt, workers)
	}

	// Step 3: Clear old currents and deposit new currents using Esirkepov algorithm
	s.Grid.ResetCurrents()
	for spIdx, sp := range s.Species {
		s.Depositor.DepositSpeciesCurrentsParallel(
			s.Grid,
			sp,
			s.oldPosBuffer[spIdx],
			dt,
			workers,
		)
	}

	// Step 4: Magnetic field second half-push: B^n -> B^{n+1/2}
	s.Solver.PushB(s.Grid, halfDt)

	// Step 5: Electric field push: E^n -> E^{n+1}
	s.Solver.PushE(s.Grid, dt)

	// Advance simulation clocks
	s.StepCount++
	s.Time += dt

	// Step 6: Diagnostic recording
	if s.Config.DiagnosticInterval > 0 && s.StepCount%s.Config.DiagnosticInterval == 0 {
		s.RecordEnergy()
	}
}

// Run executes the simulation loop for the specified number of steps.
func (s *Simulator) Run(steps int) {
	for i := 0; i < steps; i++ {
		s.Step()
	}
}
