package main

import (
	"fmt"
	"math"
	"time"

	"github.com/ZengLLQ/quasar-kinetics/pkg/diagnostics"
	"github.com/ZengLLQ/quasar-kinetics/pkg/physics"
)

func main() {
	fmt.Println("=================================================================")
	fmt.Println(" QuasarKinetics: Collisionless Landau Damping Simulation")
	fmt.Println("=================================================================")

	lc := physics.DefaultLandauConfig()
	lc.Nx = 64
	lc.NumPerCell = 80
	lc.ThermalVel = 0.25
	lc.PerturbAmp = 0.08
	lc.PerturbMode = 1

	sim, err := physics.SetupLandauSimulator(lc)
	if err != nil {
		panic(err)
	}

	k := float64(lc.PerturbMode) * (2.0 * math.Pi / lc.BoxLength)
	gammaTheory := physics.LandauTheoreticalDampingRate(k, lc.ThermalVel, 1.0)
	omegaTheory := physics.BohmGrossFrequency(k, lc.ThermalVel, 1.0)

	fmt.Printf("[+] Wavenumber k = %.4f (mode %d)\n", k, lc.PerturbMode)
	fmt.Printf("[+] Bohm-Gross frequency: omega = %.4f omega_pe\n", omegaTheory)
	fmt.Printf("[+] Theoretical Landau damping rate: gamma_L = %.4f omega_pe\n", gammaTheory)
	fmt.Println("[+] Running 150 PIC time steps...")

	sa := diagnostics.NewSpectralAnalyzer()
	initialMode := sa.DominantWaveMode(sim.Grid)
	fmt.Printf("[+] Initial Electric Field Mode 1 Amplitude: %.4e\n", initialMode.Amplitude)

	t0 := time.Now()
	sim.Run(150)
	elapsed := time.Since(t0)

	finalMode := sa.DominantWaveMode(sim.Grid)
	fmt.Printf("[+] Completed in %s\n", elapsed)
	fmt.Printf("[+] Final Electric Field Mode 1 Amplitude:   %.4e\n", finalMode.Amplitude)

	rep := diagnostics.NewEnergyReport(sim.History)
	fmt.Println("\nEnergy History:")
	fmt.Println(rep.SummaryString())
}
