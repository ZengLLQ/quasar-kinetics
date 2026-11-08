package main

import (
	"flag"
	"fmt"
	"math"
	"os"
	"time"

	"github.com/ZengLLQ/quasar-kinetics/pkg/constants"
	"github.com/ZengLLQ/quasar-kinetics/pkg/diagnostics"
	"github.com/ZengLLQ/quasar-kinetics/pkg/fields"
	"github.com/ZengLLQ/quasar-kinetics/pkg/particles"
	"github.com/ZengLLQ/quasar-kinetics/pkg/physics"
)

const banner = `
========================================================================
   ____                                   _  ___            _   _          
  / __ \__  ______ __________ _____      | |/ (_)___  ___  | |_(_)_________
 / / / / / / / __ '/ ___/ __ '/ ___/_____|   / / __ \/ _ \ | __/ / ___/ ___/
/ /_/ / /_/ / /_/ (__  ) /_/ / /  /_____/   | / / / /  __/ | |_/ / /__(__  ) 
\___\_\__,_/\__,_/____/\__,_/_/        |_|\_\/_/ /_/\___/  \__/_/\___/____/  
                                                                             
  Relativistic Particle-in-Cell (PIC) & Maxwell-Vlasov Kinetics Simulator
========================================================================
`

func main() {
	if len(os.Args) < 2 {
		printUsage()
		os.Exit(1)
	}

	subcommand := os.Args[1]

	switch subcommand {
	case "run":
		runCmd := flag.NewFlagSet("run", flag.ExitOnError)
		expFlag := runCmd.String("exp", "two-stream", "Experiment to run: two-stream, landau, cyclotron, weibel")
		stepsFlag := runCmd.Int("steps", 200, "Number of simulation time steps")
		nxFlag := runCmd.Int("nx", 64, "Number of grid cells")
		runCmd.Parse(os.Args[2:])

		executeExperiment(*expFlag, *stepsFlag, *nxFlag)

	case "benchmark":
		benchCmd := flag.NewFlagSet("benchmark", flag.ExitOnError)
		particlesFlag := benchCmd.Int("particles", 100000, "Number of particles for throughput test")
		stepsFlag := benchCmd.Int("steps", 50, "Number of pusher iterations")
		benchCmd.Parse(os.Args[2:])

		executeBenchmark(*particlesFlag, *stepsFlag)

	case "info":
		executeInfo()

	case "help", "--help", "-h":
		printUsage()

	default:
		fmt.Fprintf(os.Stderr, "Unknown subcommand: %s\n", subcommand)
		printUsage()
		os.Exit(1)
	}
}

func printUsage() {
	fmt.Print(banner)
	fmt.Println("Usage: quasar <subcommand> [flags]")
	fmt.Println("\nAvailable Subcommands:")
	fmt.Println("  run        Execute a canonical plasma kinetic experiment")
	fmt.Println("             --exp=two-stream|landau|cyclotron|weibel (default: two-stream)")
	fmt.Println("             --steps=200                             (default: 200)")
	fmt.Println("             --nx=64                                 (default: 64)")
	fmt.Println("  benchmark  Measure Boris pusher and current deposition throughput")
	fmt.Println("             --particles=100000                      (default: 100000)")
	fmt.Println("             --steps=50                              (default: 50)")
	fmt.Println("  info       Print physical constants and unit conversion systems")
	fmt.Println("  help       Display this usage instructions")
}

func executeExperiment(expName string, steps, nx int) {
	fmt.Print(banner)
	fmt.Printf("[*] Launching experiment: %s (steps=%d, nx=%d)\n", expName, steps, nx)

	switch expName {
	case "two-stream":
		tsc := physics.DefaultTwoStreamConfig()
		tsc.Nx = nx
		sim, err := physics.SetupTwoStreamSimulator(tsc)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Setup error: %v\n", err)
			os.Exit(1)
		}

		t0 := time.Now()
		fmt.Println("[+] Initialized Two-Stream Instability with counter-streaming electron beams.")
		fmt.Printf("    Theoretical maximum linear growth rate: gamma_max = %.4f omega_pe\n",
			physics.TwoStreamAnalyticalGrowthRate(1.0))

		sim.Run(steps)
		elapsed := time.Since(t0)

		fmt.Printf("\n[+] Simulation finished in %s (%d steps)\n\n", elapsed, steps)
		rep := diagnostics.NewEnergyReport(sim.History)
		fmt.Println(rep.SummaryString())

		ps := diagnostics.NewPhaseSpace2D(30, 16, 0, tsc.BoxLength, -0.4, 0.4)
		ps.BinSpecies(sim.Species[0], sim.Config.Units.SpeedOfLight)
		ps.BinSpecies(sim.Species[1], sim.Config.Units.SpeedOfLight)
		fmt.Println(ps.ASCIIHeatmap(30, 12))

	case "landau":
		lc := physics.DefaultLandauConfig()
		lc.Nx = nx
		sim, err := physics.SetupLandauSimulator(lc)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Setup error: %v\n", err)
			os.Exit(1)
		}

		k := float64(lc.PerturbMode) * (2.0 * math.Pi / lc.BoxLength)
		gamma := physics.LandauTheoreticalDampingRate(k, lc.ThermalVel, 1.0)
		omega := physics.BohmGrossFrequency(k, lc.ThermalVel, 1.0)

		fmt.Println("[+] Initialized Landau Damping simulation.")
		fmt.Printf("    Langmuir frequency: omega = %.4f\n", omega)
		fmt.Printf("    Theoretical collisionless damping rate: gamma_L = %.4f\n", gamma)

		sim.Run(steps)
		rep := diagnostics.NewEnergyReport(sim.History)
		fmt.Println(rep.SummaryString())

	case "cyclotron":
		cp := physics.CyclotronParams{
			MagneticField: 1.0,
			RelProperMom:  1.0,
			Charge:        -1.0,
			Mass:          1.0,
		}
		sim, metrics, err := physics.SetupCyclotronSimulator(cp, nx, 10.0)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Setup error: %v\n", err)
			os.Exit(1)
		}

		fmt.Println("[+] Initialized Relativistic Cyclotron orbit.")
		fmt.Printf("    Lorentz factor: gamma = %.4f\n", metrics.Gamma)
		fmt.Printf("    Relativistic gyrofrequency: omega_c = %.4f rad/s\n", metrics.OmegaC)
		fmt.Printf("    Larmor radius: r_L = %.4f\n", metrics.LarmorRadius)
		fmt.Printf("    Orbital period: T = %.4f\n", metrics.Period)

		sim.Run(steps)
		p := sim.Species[0].Particles[0]
		uMag := math.Hypot(p.Ux, p.Uy)
		fmt.Printf("[+] Finished %d steps. Proper momentum magnitude: %.6f (exact=1.000000)\n",
			steps, uMag)

	case "weibel":
		wc := physics.DefaultWeibelConfig()
		wc.Nx = nx
		sim, err := physics.SetupWeibelSimulator(wc)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Setup error: %v\n", err)
			os.Exit(1)
		}

		gamma := physics.WeibelTheoreticalGrowthRate(wc, 1.0, 1.0)
		fmt.Println("[+] Initialized Weibel Electromagnetic Filamentation Instability.")
		fmt.Printf("    Theoretical magnetic growth rate: gamma_W = %.4f omega_pe\n", gamma)

		sim.Run(steps)
		rep := diagnostics.NewEnergyReport(sim.History)
		fmt.Println(rep.SummaryString())

	default:
		fmt.Fprintf(os.Stderr, "Unknown experiment: %s\n", expName)
		os.Exit(1)
	}
}

func executeBenchmark(nParticles, steps int) {
	fmt.Print(banner)
	fmt.Printf("[*] Starting PIC Engine Benchmark with %d particles across %d steps...\n",
		nParticles, steps)

	units := constants.NewNormalizedUnits(1.0, 1.0)
	nx := 128
	dx := 0.1
	grid := fields.NewGrid1D(nx, dx, units, fields.NewPeriodicBoundary())
	pusher := particles.NewBorisPusher(units)

	sp := particles.NewSpecies(0, "bench", -1.0, 1.0, nParticles)
	for i := 0; i < nParticles; i++ {
		sp.AddParticle(particles.Particle{
			X:      float64(i%nx) * dx,
			Ux:     0.1,
			Uy:     0.05,
			Weight: 1.0,
		})
	}

	t0 := time.Now()
	for s := 0; s < steps; s++ {
		pusher.PushSpecies1DParallel(sp, grid, 0.01, 4)
	}
	elapsed := time.Since(t0)

	totalPushes := float64(nParticles) * float64(steps)
	pushesPerSec := totalPushes / elapsed.Seconds()
	megaPushes := pushesPerSec / 1e6

	fmt.Printf("[+] Benchmark Results:\n")
	fmt.Printf("    Total Particle Pushes: %.0f\n", totalPushes)
	fmt.Printf("    Elapsed Time:         %s\n", elapsed)
	fmt.Printf("    Throughput:           %.2f Million Particle-Steps / Second (Mpushes/s)\n", megaPushes)
}

func executeInfo() {
	fmt.Print(banner)
	fmt.Println("Fundamental Physical Constants (SI):")
	fmt.Printf("  Speed of light (c):             %.9e m/s\n", constants.SpeedOfLight)
	fmt.Printf("  Vacuum permittivity (eps_0):    %.9e F/m\n", constants.PermittivityOfVacuum)
	fmt.Printf("  Vacuum permeability (mu_0):     %.9e N/A^2\n", constants.PermeabilityOfVacuum)
	fmt.Printf("  Elementary charge (e):          %.9e C\n", constants.ElementaryCharge)
	fmt.Printf("  Electron rest mass (m_e):       %.9e kg\n", constants.ElectronMass)
	fmt.Printf("  Proton rest mass (m_p):         %.9e kg\n", constants.ProtonMass)
	fmt.Printf("  Boltzmann constant (k_B):       %.9e J/K\n", constants.BoltzmannConstant)
	fmt.Println("\nNormalized Unit System Scaling:")
	fmt.Println("  c = 1, epsilon_0 = 1, mu_0 = 1, m_e = 1, q_e = 1")
	fmt.Println("  Length scale: skin depth c / omega_pe")
	fmt.Println("  Time scale:   inverse plasma frequency 1 / omega_pe")
}
