package main

import (
	"fmt"
	"time"

	"github.com/ZengLLQ/quasar-kinetics/pkg/diagnostics"
	"github.com/ZengLLQ/quasar-kinetics/pkg/physics"
)

func main() {
	fmt.Println("=================================================================")
	fmt.Println(" QuasarKinetics: Two-Stream Instability Simulation")
	fmt.Println("=================================================================")

	tsc := physics.DefaultTwoStreamConfig()
	tsc.Nx = 64
	tsc.NumPerCell = 40
	tsc.PerturbAmp = 0.02
	tsc.DriftVel = 0.2

	sim, err := physics.SetupTwoStreamSimulator(tsc)
	if err != nil {
		panic(err)
	}

	gammaTheory := physics.TwoStreamAnalyticalGrowthRate(1.0)
	fmt.Printf("[+] Drift velocity: v0 = +/-%.2f c\n", tsc.DriftVel)
	fmt.Printf("[+] Theoretical maximum growth rate: gamma_max = %.4f omega_pe\n", gammaTheory)
	fmt.Println("[+] Running 200 PIC time steps...")

	t0 := time.Now()
	sim.Run(200)
	elapsed := time.Since(t0)

	fmt.Printf("[+] Completed in %s\n\n", elapsed)

	rep := diagnostics.NewEnergyReport(sim.History)
	fmt.Println("Energy Evolution:")
	fmt.Println(rep.SummaryString())

	fmt.Println("Final Phase-Space Vortex Distribution f(x, vx):")
	ps := diagnostics.NewPhaseSpace2D(32, 16, 0, tsc.BoxLength, -0.45, 0.45)
	ps.BinSpecies(sim.Species[0], sim.Config.Units.SpeedOfLight)
	ps.BinSpecies(sim.Species[1], sim.Config.Units.SpeedOfLight)
	fmt.Println(ps.ASCIIHeatmap(32, 12))
}
