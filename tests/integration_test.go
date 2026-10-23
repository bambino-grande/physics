package tests

import (
	"math"
	"nbody/physics"
	"testing"
)

func TestIntegrationFullRun(t *testing.T) {
	const (
		n     = 600
		steps = 400
		dt    = 0.002
	)
	for _, name := range []string{physics.Uniform, physics.DenseCore, physics.Collision} {
		t.Run(name, func(t *testing.T) {
			b, err := physics.Generate(name, n, 11, tG, tEps)
			if err != nil {
				t.Fatal(err)
			}
			s := physics.BarnesHut{G: tG, Eps: tEps, Theta: 0.5, Workers: 4}
			e0 := physics.TotalEnergy(b, tG, tEps, 4)
			physics.Init(b, s)
			var inter int64
			for i := 0; i < steps; i++ {
				inter += physics.Step(b, dt, s).Interactions
			}
			if len(b) != n {
				t.Fatalf("число тел изменилось: %d", len(b))
			}
			var mass float64
			for _, p := range b {
				mass += p.Mass
				if math.IsNaN(p.Pos.X+p.Pos.Y+p.Vel.X+p.Vel.Y) || math.IsInf(p.Pos.X+p.Pos.Y+p.Vel.X+p.Vel.Y, 0) {
					t.Fatalf("тело %d: NaN/Inf", p.ID)
				}
			}
			if math.Abs(mass-1) > 1e-12 {
				t.Errorf("суммарная масса %.15f", mass)
			}
			de := math.Abs(physics.TotalEnergy(b, tG, tEps, 4)-e0) / math.Abs(e0)
			p := physics.Momentum(b).Norm()
			t.Logf("%s: |ΔE/E0|=%.2e |P|=%.2e взаимодействий/шаг=%d", name, de, p, inter/steps)
			if de > 0.05 {
				t.Errorf("дрейф энергии %.2e", de)
			}
			if p > 0.05 {
				t.Errorf("импульс ушёл: %.2e", p)
			}
			if inter == 0 {
				t.Error("взаимодействия не считались")
			}
		})
	}
}

func TestIntegrationBHMatchesDirectTrajectory(t *testing.T) {
	b, _ := physics.Generate(physics.DenseCore, 400, 12, tG, tEps)
	d, h := clone(b), clone(b)
	sd := physics.Direct{G: tG, Eps: tEps, Workers: 2}
	sh := physics.BarnesHut{G: tG, Eps: tEps, Theta: 0.3, Workers: 2}
	physics.Init(d, sd)
	physics.Init(h, sh)
	for i := 0; i < 100; i++ {
		physics.Step(d, 0.002, sd)
		physics.Step(h, 0.002, sh)
	}
	var worst float64
	for i := range d {
		worst = math.Max(worst, d[i].Pos.Sub(h[i].Pos).Norm())
	}
	t.Logf("max расхождение позиций после 100 шагов: %.2e", worst)

	if worst > 0.02 {
		t.Errorf("траектории BH и прямого расчёта разошлись: %.2e", worst)
	}
}
