package tests

import (
	"math"
	"nbody/physics"
	"nbody/physics/dto"
	"sort"
	"testing"
)

const (
	tG   = 1.0
	tEps = 0.02
)

func relErrs(a, ref []dto.Body) []float64 {
	e := make([]float64, len(a))
	for i := range a {
		e[i] = a[i].Acc.Sub(ref[i].Acc).Norm() / ref[i].Acc.Norm()
	}
	sort.Float64s(e)
	return e
}

func clone(b []dto.Body) []dto.Body { return append([]dto.Body(nil), b...) }

func TestTwoBodyCircularOrbit(t *testing.T) {

	m, d := 1.0, 1.0
	v := math.Sqrt(tG * m / (2 * d))
	T := 2 * math.Pi * (d / 2) / v
	b := []dto.Body{
		{ID: 0, Mass: m, Pos: dto.Vec2{X: -d / 2, Y: 0}, Vel: dto.Vec2{X: 0, Y: -v}},
		{ID: 1, Mass: m, Pos: dto.Vec2{X: d / 2, Y: 0}, Vel: dto.Vec2{X: 0, Y: v}},
	}
	start := clone(b)
	s := physics.Direct{G: tG, Eps: 0}
	e0 := physics.TotalEnergy(b, tG, 0, 1)
	steps := 4000
	dt := T / float64(steps)
	physics.Init(b, s)
	for i := 0; i < steps; i++ {
		physics.Step(b, dt, s)
	}
	for i := range b {
		if dist := b[i].Pos.Sub(start[i].Pos).Norm(); dist > 1e-4 {
			t.Errorf("тело %d не вернулось за период: отклонение %.2e", i, dist)
		}
	}
	if de := math.Abs(physics.TotalEnergy(b, tG, 0, 1)-e0) / math.Abs(e0); de > 1e-6 {
		t.Errorf("дрейф энергии %.2e", de)
	}
}

func TestMomentumConservedDirect(t *testing.T) {
	b, _ := physics.Generate(physics.DenseCore, 200, 1, tG, tEps)
	s := physics.Direct{G: tG, Eps: tEps}
	physics.Init(b, s)
	for i := 0; i < 200; i++ {
		physics.Step(b, 0.002, s)
	}
	if p := physics.Momentum(b).Norm(); p > 1e-12 {
		t.Errorf("импульс не сохраняется: |P| = %.2e", p)
	}
}

func TestBarnesHutThetaZeroEqualsDirect(t *testing.T) {
	b, _ := physics.Generate(physics.DenseCore, 500, 2, tG, tEps)
	ref, bh := clone(b), clone(b)
	physics.Direct{G: tG, Eps: tEps}.Accel(ref)
	physics.BarnesHut{G: tG, Eps: tEps, Theta: 0}.Accel(bh)
	e := relErrs(bh, ref)
	if e[len(e)-1] > 1e-10 {
		t.Errorf("theta=0 должна совпадать с прямым расчётом, max ошибка %.2e", e[len(e)-1])
	}
}

func TestBarnesHutErrorGrowsWithTheta(t *testing.T) {
	b, _ := physics.Generate(physics.DenseCore, 3000, 3, tG, tEps)
	ref := clone(b)
	physics.Direct{G: tG, Eps: tEps, Workers: 2}.Accel(ref)
	prev := 0.0
	for _, th := range []float64{0.2, 0.5, 1.0} {
		x := clone(b)
		physics.BarnesHut{G: tG, Eps: tEps, Theta: th}.Accel(x)
		e := relErrs(x, ref)
		med := e[len(e)/2]
		t.Logf("theta=%.1f  медианная ошибка силы %.2e", th, med)
		if med <= prev {
			t.Errorf("ошибка должна расти с theta: %.2e после %.2e", med, prev)
		}
		prev = med
		if th == 0.5 && med > 0.01 {
			t.Errorf("при theta=0.5 медианная ошибка должна быть < 1%%, получили %.2e", med)
		}
	}
}

func TestBarnesHutFewerInteractions(t *testing.T) {
	b, _ := physics.Generate(physics.Uniform, 4000, 4, tG, tEps)
	d := physics.Direct{G: tG, Eps: tEps}.Accel(clone(b))
	h := physics.BarnesHut{G: tG, Eps: tEps, Theta: 0.5}.Accel(clone(b))
	if h.Interactions >= d.Interactions/5 {
		t.Errorf("BH должен считать заметно меньше взаимодействий: %d против %d", h.Interactions, d.Interactions)
	}
}

func TestParallelEqualsSerial(t *testing.T) {
	b, _ := physics.Generate(physics.Collision, 1000, 5, tG, tEps)
	s1, s4 := clone(b), clone(b)
	physics.BarnesHut{G: tG, Eps: tEps, Theta: 0.5, Workers: 1}.Accel(s1)
	physics.BarnesHut{G: tG, Eps: tEps, Theta: 0.5, Workers: 4}.Accel(s4)
	for i := range s1 {
		if s1[i].Acc != s4[i].Acc {
			t.Fatalf("тело %d: многопоточный результат отличается от однопоточного", i)
		}
	}
}

func TestLeapfrogTimeReversible(t *testing.T) {
	b, _ := physics.Generate(physics.DenseCore, 150, 6, tG, tEps)
	start := clone(b)
	s := physics.Direct{G: tG, Eps: tEps}
	physics.Init(b, s)
	for i := 0; i < 100; i++ {
		physics.Step(b, 0.002, s)
	}
	for i := range b {
		b[i].Vel = b[i].Vel.Scale(-1)
	}
	physics.Init(b, s)
	for i := 0; i < 100; i++ {
		physics.Step(b, 0.002, s)
	}
	var worst float64
	for i := range b {
		worst = math.Max(worst, b[i].Pos.Sub(start[i].Pos).Norm())
	}
	if worst > 1e-8 {
		t.Errorf("обратимость нарушена: max отклонение %.2e", worst)
	}
}

func TestEnergyDrift(t *testing.T) {
	b, _ := physics.Generate(physics.DenseCore, 300, 7, tG, tEps)
	s := physics.Direct{G: tG, Eps: tEps}
	e0 := physics.TotalEnergy(b, tG, tEps, 1)
	physics.Init(b, s)
	for i := 0; i < 500; i++ {
		physics.Step(b, 0.001, s)
	}
	de := math.Abs(physics.TotalEnergy(b, tG, tEps, 1)-e0) / math.Abs(e0)
	t.Logf("|ΔE/E0| за 500 шагов = %.2e", de)
	if de > 1e-3 {
		t.Errorf("слишком большой дрейф энергии: %.2e", de)
	}
}

func TestTreeMassAndCenter(t *testing.T) {
	b, _ := physics.Generate(physics.Collision, 800, 8, tG, tEps)
	tr := physics.BuildTree(b)
	if math.Abs(tr.RootMass()-1) > 1e-12 {
		t.Errorf("масса корня %.15f, ожидалась 1", tr.RootMass())
	}
	if c := tr.RootCOM().Norm(); c > 1e-12 {
		t.Errorf("центр масс должен быть в нуле, получили %.2e", c)
	}
}

func TestCoincidentBodiesDoNotHang(t *testing.T) {
	b := []dto.Body{
		{ID: 0, Mass: 1, Pos: dto.Vec2{X: 0.3, Y: 0.3}},
		{ID: 1, Mass: 1, Pos: dto.Vec2{X: 0.3, Y: 0.3}},
		{ID: 2, Mass: 1, Pos: dto.Vec2{X: -1, Y: 2}},
	}
	physics.BarnesHut{G: tG, Eps: tEps, Theta: 0.5}.Accel(b)
	for _, x := range b {
		if math.IsNaN(x.Acc.X) || math.IsInf(x.Acc.X, 0) {
			t.Fatal("NaN/Inf в ускорении")
		}
	}
}

func TestScenariosDeterministicAndZeroDrift(t *testing.T) {
	for _, name := range []string{physics.Uniform, physics.DenseCore, physics.Collision} {
		a, _ := physics.Generate(name, 500, 42, tG, tEps)
		c, _ := physics.Generate(name, 500, 42, tG, tEps)
		for i := range a {
			if a[i] != c[i] {
				t.Fatalf("%s: генерация не детерминирована", name)
			}
		}
		if p := physics.Momentum(a).Norm(); p > 1e-12 {
			t.Errorf("%s: суммарный импульс %.2e", name, p)
		}
	}
	if _, err := physics.Generate("nope", 10, 1, tG, tEps); err == nil {
		t.Error("ожидалась ошибка для неизвестного сценария")
	}
}
