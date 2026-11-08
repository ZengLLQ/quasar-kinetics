package diagnostics

import (
	"fmt"
	"strings"

	"github.com/ZengLLQ/quasar-kinetics/pkg/engine"
)

// EnergyReport formats simulation energy history for display and CSV serialization.
type EnergyReport struct {
	Snapshots []engine.EnergySnapshot
}

// NewEnergyReport wraps the snapshot history.
func NewEnergyReport(history []engine.EnergySnapshot) *EnergyReport {
	return &EnergyReport{Snapshots: history}
}

// ToCSV serializes energy history to comma-separated values format.
func (er *EnergyReport) ToCSV() string {
	var sb strings.Builder
	sb.WriteString("step,time,kinetic_energy,electric_energy,magnetic_energy,total_energy,rel_drift\n")
	for _, s := range er.Snapshots {
		sb.WriteString(fmt.Sprintf("%d,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
			s.Step, s.Time, s.KineticEnergy, s.ElectricEnergy, s.MagneticEnergy, s.TotalEnergy, s.RelDrift))
	}
	return sb.String()
}

// SummaryString creates a compact human-readable ASCII summary table.
func (er *EnergyReport) SummaryString() string {
	var sb strings.Builder
	sb.WriteString(fmt.Sprintf("%-6s | %-10s | %-12s | %-12s | %-12s | %-12s | %-10s\n",
		"Step", "Time", "E_kin", "E_elec", "E_mag", "E_total", "Drift"))
	sb.WriteString(strings.Repeat("-", 85) + "\n")

	// Print up to 10 spaced samples
	n := len(er.Snapshots)
	if n == 0 {
		return "No energy records.\n"
	}

	stride := 1
	if n > 10 {
		stride = n / 10
	}

	for i := 0; i < n; i += stride {
		s := er.Snapshots[i]
		sb.WriteString(fmt.Sprintf("%-6d | %-10.4f | %-12.4e | %-12.4e | %-12.4e | %-12.4e | %-10.2e\n",
			s.Step, s.Time, s.KineticEnergy, s.ElectricEnergy, s.MagneticEnergy, s.TotalEnergy, s.RelDrift))
	}

	// Always print the final step if not already printed
	lastIdx := n - 1
	if (lastIdx % stride) != 0 {
		s := er.Snapshots[lastIdx]
		sb.WriteString(fmt.Sprintf("%-6d | %-10.4f | %-12.4e | %-12.4e | %-12.4e | %-12.4e | %-10.2e\n",
			s.Step, s.Time, s.KineticEnergy, s.ElectricEnergy, s.MagneticEnergy, s.TotalEnergy, s.RelDrift))
	}

	return sb.String()
}
