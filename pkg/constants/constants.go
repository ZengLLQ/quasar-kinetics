// Package constants defines fundamental physical constants in SI units
// and provides unit systems for relativistic particle-in-cell simulations.
package constants

import "math"

// Fundamental physical constants in SI units.
const (
	// SpeedOfLight is the speed of light in vacuum c [m/s].
	SpeedOfLight = 299792458.0

	// PermittivityOfVacuum is the vacuum electric permittivity epsilon_0 [F/m].
	PermittivityOfVacuum = 8.8541878128e-12

	// PermeabilityOfVacuum is the vacuum magnetic permeability mu_0 [N/A^2].
	PermeabilityOfVacuum = 1.25663706212e-6

	// ElementaryCharge is the proton charge magnitude e [C].
	ElementaryCharge = 1.602176634e-19

	// ElectronMass is the electron rest mass m_e [kg].
	ElectronMass = 9.1093837015e-31

	// ProtonMass is the proton rest mass m_p [kg].
	ProtonMass = 1.67262192369e-27

	// BoltzmannConstant is the Boltzmann constant k_B [J/K].
	BoltzmannConstant = 1.380649e-23
)

// UnitSystem defines the dimensional or dimensionless scaling for PIC simulations.
type UnitSystem int

const (
	// UnitSI uses standard International System of Units (m, s, kg, C, V, T).
	UnitSI UnitSystem = iota
	// UnitNormalized uses dimensionless plasma units where c = 1, eps_0 = 1, mu_0 = 1, m_e = 1, q_e = 1.
	UnitNormalized
)

// Units holds scaling parameters for converting between physical SI and internal normalized units.
type Units struct {
	System         UnitSystem
	SpeedOfLight   float64
	Permittivity   float64
	Permeability   float64
	ElementaryMass float64
	ElementaryQ    float64

	// Characteristic scales for conversion to SI
	LengthScale   float64 // [m] per simulation unit
	TimeScale     float64 // [s] per simulation unit
	DensityScale  float64 // [m^-3] per simulation unit
	MassScale     float64 // [kg] per simulation unit
	ChargeScale   float64 // [C] per simulation unit
}

// NewSIUnits returns a Units struct configured with standard SI physical constants.
func NewSIUnits() *Units {
	return &Units{
		System:         UnitSI,
		SpeedOfLight:   SpeedOfLight,
		Permittivity:   PermittivityOfVacuum,
		Permeability:   PermeabilityOfVacuum,
		ElementaryMass: ElectronMass,
		ElementaryQ:    ElementaryCharge,
		LengthScale:    1.0,
		TimeScale:      1.0,
		DensityScale:   1.0,
		MassScale:      1.0,
		ChargeScale:    1.0,
	}
}

// NewNormalizedUnits returns dimensionless units with c=1, eps_0=1, mu_0=1, m_e=1, q_e=1.
// plasmaOmega is the reference electron plasma frequency omega_pe [rad/s].
// referenceDensity is the reference number density n_0 [m^-3].
func NewNormalizedUnits(plasmaOmega, referenceDensity float64) *Units {
	if plasmaOmega <= 0 {
		plasmaOmega = 1.0e12 // Default 1 THz
	}
	if referenceDensity <= 0 {
		referenceDensity = 1.0e20 // Default 10^20 m^-3
	}

	timeScale := 1.0 / plasmaOmega
	lengthScale := SpeedOfLight / plasmaOmega // Skin depth c / omega_pe
	chargeScale := ElementaryCharge
	massScale := ElectronMass

	return &Units{
		System:         UnitNormalized,
		SpeedOfLight:   1.0,
		Permittivity:   1.0,
		Permeability:   1.0,
		ElementaryMass: 1.0,
		ElementaryQ:    1.0,
		LengthScale:    lengthScale,
		TimeScale:      timeScale,
		DensityScale:   referenceDensity,
		MassScale:      massScale,
		ChargeScale:    chargeScale,
	}
}

// PlasmaFrequency calculates the cold electron plasma frequency:
//
//	omega_pe = sqrt( n * e^2 / (epsilon_0 * m_e) )
func PlasmaFrequency(n float64, u *Units) float64 {
	return math.Sqrt(n * math.Pow(u.ElementaryQ, 2) / (u.Permittivity * u.ElementaryMass))
}

// DebyeLength calculates the electron Debye shielding length:
//
//	lambda_D = sqrt( epsilon_0 * k_B * T / (n * e^2) )
func DebyeLength(n, temperatureK float64, u *Units) float64 {
	kbT := BoltzmannConstant * temperatureK
	if u.System == UnitNormalized {
		// In normalized units, temperature is provided directly in thermal energy units
		kbT = temperatureK
	}
	return math.Sqrt((u.Permittivity * kbT) / (n * math.Pow(u.ElementaryQ, 2)))
}

// SkinDepth calculates the collisionless plasma skin depth:
//
//	d_e = c / omega_pe
func SkinDepth(n float64, u *Units) float64 {
	wpe := PlasmaFrequency(n, u)
	if wpe == 0 {
		return math.Inf(1)
	}
	return u.SpeedOfLight / wpe
}

// RelativisticGamma computes the Lorentz factor:
//
//	gamma = 1 / sqrt(1 - (v/c)^2) = sqrt(1 + (u/c)^2)
//
// where u = gamma * v is the normalized 4-momentum / proper velocity.
func RelativisticGamma(ux, uy, uz, c float64) float64 {
	uSq := (ux*ux + uy*uy + uz*uz) / (c * c)
	return math.Sqrt(1.0 + uSq)
}

// VelocityFromProperMom computes spatial velocity v = u / gamma.
func VelocityFromProperMom(ux, uy, uz, c float64) (vx, vy, vz float64) {
	gamma := RelativisticGamma(ux, uy, uz, c)
	return ux / gamma, uy / gamma, uz / gamma
}
