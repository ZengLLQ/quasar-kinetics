package constants_test

import (
	"math"
	"testing"

	"github.com/ZengLLQ/quasar-kinetics/pkg/constants"
)

func TestConstantsElectromagneticRelations(t *testing.T) {
	// In SI: c = 1 / sqrt(epsilon_0 * mu_0)
	cDerived := 1.0 / math.Sqrt(constants.PermittivityOfVacuum*constants.PermeabilityOfVacuum)
	relErr := math.Abs(cDerived-constants.SpeedOfLight) / constants.SpeedOfLight
	if relErr > 1e-4 {
		t.Errorf("SI electromagnetic relation mismatch: derived c=%g, expected %g, relErr=%g",
			cDerived, constants.SpeedOfLight, relErr)
	}
}

func TestUnitsCreation(t *testing.T) {
	si := constants.NewSIUnits()
	if si.SpeedOfLight != constants.SpeedOfLight {
		t.Fatalf("expected SI c to be %g, got %g", constants.SpeedOfLight, si.SpeedOfLight)
	}
	if si.System != constants.UnitSI {
		t.Fatalf("expected UnitSI, got %v", si.System)
	}

	norm := constants.NewNormalizedUnits(1e12, 1e20)
	if norm.SpeedOfLight != 1.0 || norm.Permittivity != 1.0 {
		t.Fatalf("expected normalized c=1 and eps0=1, got c=%g, eps0=%g",
			norm.SpeedOfLight, norm.Permittivity)
	}
}

func TestPlasmaFrequencyAndDebyeLength(t *testing.T) {
	norm := constants.NewNormalizedUnits(1.0, 1.0)
	// In normalized units with n=1, e=1, eps=1, m=1 => omega_pe = 1
	wp := constants.PlasmaFrequency(1.0, norm)
	if math.Abs(wp-1.0) > 1e-12 {
		t.Errorf("expected omega_pe = 1.0, got %g", wp)
	}

	// Debye length with kbT = 0.01, n = 1, e = 1 => lambda_D = sqrt(0.01) = 0.1
	ld := constants.DebyeLength(1.0, 0.01, norm)
	if math.Abs(ld-0.1) > 1e-12 {
		t.Errorf("expected lambda_D = 0.1, got %g", ld)
	}

	// Skin depth with c=1, omega_pe=1 => d_e = 1.0
	sd := constants.SkinDepth(1.0, norm)
	if math.Abs(sd-1.0) > 1e-12 {
		t.Errorf("expected skin depth = 1.0, got %g", sd)
	}
}

func TestRelativisticGamma(t *testing.T) {
	c := 1.0
	// Non-relativistic: u = 0 => gamma = 1
	g0 := constants.RelativisticGamma(0, 0, 0, c)
	if g0 != 1.0 {
		t.Errorf("expected gamma=1 for stationary particle, got %g", g0)
	}

	// Relativistic: u = c => gamma = sqrt(1 + 1) = sqrt(2) approx 1.41421356
	g1 := constants.RelativisticGamma(1.0, 0, 0, c)
	if math.Abs(g1-math.Sqrt(2)) > 1e-12 {
		t.Errorf("expected gamma=sqrt(2), got %g", g1)
	}

	// Check velocity inversion: v = u / gamma => v < c always
	vx, vy, vz := constants.VelocityFromProperMom(10.0, 0, 0, c)
	if vx >= c {
		t.Errorf("velocity must be strictly subluminal, got vx=%g (c=%g)", vx, c)
	}
	if vy != 0 || vz != 0 {
		t.Errorf("transverse velocities must be zero, got (%g, %g)", vy, vz)
	}
}
