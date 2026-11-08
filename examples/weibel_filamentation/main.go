package main

import (
	"fmt"
	"time"

	"github.com/ZengLLQ/quasar-kinetics/pkg/diagnostics"
	"github.com/ZengLLQ/quasar-kinetics/pkg/physics"
)

func main() {
	fmt.Println("=================================================================")
	fmt.Println(" QuasarKinetics: Weibel Electromagnetic Filamentation")
	fmt.Println("=================================================================")

	wc := physics.DefaultWeibelConfig()
	wc.Nx = 64
	wc.NumPerCell = 60
	wc.VthX = 0.04
	wc.VthY = 0.28 // Strong thermal anisotropy Ty / Tx approx 49

	sim, err := physics.SetupWeibelSimulator(wc)
	if err != nil {
		panic(err)
	}

	gammaW := physics.WeibelTheoreticalGrowthRate(wc, 1.0, 1.0)
	fmt.Printf("[+] Temperature Anisotropy: (Ty - Tx) / Ty = %.2f%%\n",
		100.0*(wc.VthY*wc.VthY-wc.VthX*wc.VthX)/(wc.VthY*wc.VthY))
	fmt.Printf("[+] Theoretical Weibel Growth Rate: gamma_W = %.4f omega_pe\n", gammaW)
	fmt.Printf("[+] Initial Magnetic Energy:        %.4e\n", sim.Grid.MagneticEnergy())

	t0 := time.Now()
	sim.Run(200)
	elapsed := time.Since(t0)

	fmt.Printf("[+] Completed in %s\n", elapsed)
	fmt.Printf("[+] Final Magnetic Energy:          %.4e\n", sim.Grid.MagneticEnergy())

	rep := diagnostics.NewEnergyReport(sim.History)
	fmt.Println("\nEnergy Evolution:")
	fmt.Println(rep.SummaryString())
}
