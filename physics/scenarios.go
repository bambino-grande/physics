package physics

import (
	"fmt"
	"math"
	"math/rand"
	"nbody/physics/dto"
	"sort"
)

const (
	Uniform   = "uniform"
	DenseCore = "dense_core"
	Collision = "collision"
)

func Generate(name string, n int, seed int64, g, eps float64) ([]dto.Body, error) {
	rng := rand.New(rand.NewSource(seed))
	var b []dto.Body
	switch name {
	case Uniform:
		b = make([]dto.Body, n)
		for i := range b {
			r := math.Sqrt(rng.Float64())
			a := 2 * math.Pi * rng.Float64()
			b[i].Pos = dto.Vec2{X: r * math.Cos(a), Y: r * math.Sin(a)}
		}
		circular(b, dto.Vec2{}, g, eps, 1)
	case DenseCore:
		b = make([]dto.Body, n)
		nCore := n * 8 / 10
		for i := range b {
			if i < nCore {
				b[i].Pos = dto.Vec2{X: rng.NormFloat64() * 0.05, Y: rng.NormFloat64() * 0.05}
			} else {
				r := math.Sqrt(rng.Float64())
				a := 2 * math.Pi * rng.Float64()
				b[i].Pos = dto.Vec2{X: r * math.Cos(a), Y: r * math.Sin(a)}
			}
		}
		circular(b, dto.Vec2{}, g, eps, 1)
	case Collision:
		half := n / 2
		b = make([]dto.Body, n)
		for i := range b {
			b[i].Pos = dto.Vec2{X: rng.NormFloat64() * 0.1, Y: rng.NormFloat64() * 0.1}
		}
		circular(b[:half], dto.Vec2{}, g, eps, 0.5)
		circular(b[half:], dto.Vec2{}, g, eps, 0.5)
		shift(b[:half], dto.Vec2{X: -0.6, Y: 0.05}, dto.Vec2{X: 0.5, Y: 0})
		shift(b[half:], dto.Vec2{X: 0.6, Y: -0.05}, dto.Vec2{X: -0.5, Y: 0})
	default:
		return nil, fmt.Errorf("неизвестный сценарий %q (есть: %s, %s, %s)", name, Uniform, DenseCore, Collision)
	}
	m := 1.0 / float64(n)
	for i := range b {
		b[i].ID = i
		b[i].Mass = m
	}
	removeDrift(b)
	return b, nil
}

func circular(b []dto.Body, center dto.Vec2, g, eps, totalMass float64) {
	n := len(b)
	idx := make([]int, n)
	for i := range idx {
		idx[i] = i
	}
	rad := func(i int) float64 { return b[i].Pos.Sub(center).Norm() }
	sort.Slice(idx, func(x, y int) bool { return rad(idx[x]) < rad(idx[y]) })
	rng := rand.New(rand.NewSource(int64(n)*7919 + 1))
	for k, i := range idx {
		d := b[i].Pos.Sub(center)
		r := d.Norm()
		if r == 0 {
			continue
		}
		mEnc := totalMass * float64(k+1) / float64(n)
		v := math.Sqrt(g * mEnc * r * r / math.Pow(r*r+eps*eps, 1.5))
		v *= 1 + 0.1*rng.NormFloat64()
		b[i].Vel = dto.Vec2{X: -d.Y / r * v, Y: d.X / r * v}
	}
}

func shift(b []dto.Body, dp, dv dto.Vec2) {
	for i := range b {
		b[i].Pos = b[i].Pos.Add(dp)
		b[i].Vel = b[i].Vel.Add(dv)
	}
}

func removeDrift(b []dto.Body) {
	var cm, p dto.Vec2
	var m float64
	for _, x := range b {
		cm = cm.Add(x.Pos.Scale(x.Mass))
		p = p.Add(x.Vel.Scale(x.Mass))
		m += x.Mass
	}
	cm, p = cm.Scale(1/m), p.Scale(1/m)
	for i := range b {
		b[i].Pos = b[i].Pos.Sub(cm)
		b[i].Vel = b[i].Vel.Sub(p)
	}
}
