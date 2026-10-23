package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"time"

	"nbody/physics"
)

func main() {
	scenario := flag.String("scenario", physics.DenseCore, "uniform | dense_core | collision")
	n := flag.Int("n", 2000, "число тел")
	steps := flag.Int("steps", 1000, "число шагов")
	dt := flag.Float64("dt", 0.002, "шаг по времени")
	eps := flag.Float64("eps", 0.02, "сглаживание ε")
	theta := flag.Float64("theta", 0.5, "параметр θ Барнса—Хата")
	solver := flag.String("solver", "bh", "bh | direct")
	workers := flag.Int("workers", physics.NumCPU(), "число потоков")
	seed := flag.Int64("seed", 1, "seed генератора")
	every := flag.Int("every", 100, "печатать диагностику каждые N шагов")
	energy := flag.Bool("energy", true, "считать полную энергию (O(N²), только в точках вывода)")
	logPath := flag.String("log", "", "CSV с диагностикой по шагам вывода")
	snapPath := flag.String("snap", "", "CSV со снимками позиций (step,id,x,y) для визуализации")
	snapEvery := flag.Int("snap-every", 50, "как часто писать снимок")
	flag.Parse()

	const G = 1.0
	b, err := physics.Generate(*scenario, *n, *seed, G, *eps)
	if err != nil {
		log.Fatal(err)
	}
	var s physics.Solver
	switch *solver {
	case "bh":
		s = physics.BarnesHut{G: G, Eps: *eps, Theta: *theta, Workers: *workers}
	case "direct":
		s = physics.Direct{G: G, Eps: *eps, Workers: *workers}
	default:
		log.Fatalf("неизвестный solver %q", *solver)
	}

	var logW, snapW *bufio.Writer
	if *logPath != "" {
		f, err := os.Create(*logPath)
		if err != nil {
			log.Fatal(err)
		}
		defer f.Close()
		logW = bufio.NewWriter(f)
		defer logW.Flush()
		fmt.Fprintln(logW, "step,t,energy,dE_over_E0,momentum,angular_momentum,interactions,step_ms")
	}
	if *snapPath != "" {
		f, err := os.Create(*snapPath)
		if err != nil {
			log.Fatal(err)
		}
		defer f.Close()
		snapW = bufio.NewWriter(f)
		defer snapW.Flush()
		fmt.Fprintln(snapW, "step,id,x,y")
	}

	var e0 float64
	if *energy {
		e0 = physics.TotalEnergy(b, G, *eps, *workers)
	}
	physics.Init(b, s)
	fmt.Printf("сценарий=%s N=%d solver=%s theta=%.2f dt=%g eps=%g workers=%d seed=%d\n", *scenario, *n, *solver, *theta, *dt, *eps, *workers, *seed)
	fmt.Printf("%7s %9s %13s %11s %10s %12s %12s %9s\n", "step", "t", "E", "dE/E0", "|P|", "Lz", "interactions", "step_ms")

	var stepTime time.Duration
	var inter int64
	for k := 1; k <= *steps; k++ {
		t0 := time.Now()
		st := physics.Step(b, *dt, s)
		stepTime += time.Since(t0)
		inter += st.Interactions
		if snapW != nil && k%*snapEvery == 0 {
			for _, p := range b {
				fmt.Fprintf(snapW, "%d,%d,%.6f,%.6f\n", k, p.ID, p.Pos.X, p.Pos.Y)
			}
		}
		if k%*every == 0 || k == *steps {
			ms := float64(stepTime.Microseconds()) / 1000 / float64(*every)
			var e, de float64
			if *energy {
				e = physics.TotalEnergy(b, G, *eps, *workers)
				de = math.Abs(e-e0) / math.Abs(e0)
			}
			p, lz := physics.Momentum(b).Norm(), physics.AngularMomentum(b)
			fmt.Printf("%7d %9.4f %13.6f %11.2e %10.2e %12.6f %12d %9.3f\n", k, float64(k)**dt, e, de, p, lz, inter/int64(*every), ms)
			if logW != nil {
				fmt.Fprintf(logW, "%d,%.6f,%.9f,%.3e,%.3e,%.9f,%d,%.4f\n", k, float64(k)**dt, e, de, p, lz, inter/int64(*every), ms)
			}
			stepTime, inter = 0, 0
		}
	}
}
