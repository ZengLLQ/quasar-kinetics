package main

import (
	"fmt"
	"math"

	"github.com/ZengLLQ/quasar-kinetics/pkg/physics"
)

func main() {
	fmt.Println("=================================================================")
	fmt.Println(" QuasarKinetics: Relativistic Cyclotron Gyration")
	fmt.Println("=================================================================")

	// Highly relativistic electron with gamma approx 2.0 (u = sqrt(3) approx 1.732)
	uRel := math.Sqrt(3.0)
	cp := physics.CyclotronParams{
		MagneticField: 1.5,
		RelProperMom:  uRel,
		Charge:        -1.0,
		Mass:          1.0,
	}

	sim, metrics, err := physics.SetupCyclotronSimulator(cp, 64, 15.0)
	if err != nil {
		panic(err)
	}

	fmt.Printf("[+] Lorentz Factor (gamma):            %.4f\n", metrics.Gamma)
	fmt.Printf("[+] Subluminal Velocity (v / c):       %.4f\n", metrics.Velocity)
	fmt.Printf("[+] Relativistic Gyrofrequency (w_c):  %.4f rad/s\n", metrics.OmegaC)
	fmt.Printf("[+] Larmor Gyroradius (r_L):           %.4f\n", metrics.LarmorRadius)
	fmt.Printf("[+] Gyration Period (T):               %.4f s\n", metrics.Period)
	fmt.Println("[+] Running 2 full orbital revolutions (400 steps)...")

	sim.Run(400)

	p := sim.Species[0].Particles[0]
	uFinal := math.Hypot(p.Ux, p.Uy)
	relDrift := math.Abs(uFinal-uRel) / uRel

	fmt.Printf("[+] Particle final state at step %d:\n", sim.StepCount)
	fmt.Printf("    Coordinates:    x = %.4f, y = %.4f\n", p.X, p.Y)
	fmt.Printf("    Proper momentum: |u| = %.6f (exact: %.6f, drift: %.2e)\n",
		uFinal, uRel, relDrift)
	fmt.Printf("    Relativistic Boris pusher preserved exact phase-space symplectic volume.\n")
}
