package physics

import (
	"math"
	"nbody/physics/dto"
)

func Kinetic(b []dto.Body) float64 {
	var e float64
	for _, p := range b {
		e += 0.5 * p.Mass * (p.Vel.X*p.Vel.X + p.Vel.Y*p.Vel.Y)
	}
	return e
}

func Potential(b []dto.Body, g, eps float64, workers int) float64 {
	eps2 := eps * eps
	n := len(b)
	partial := make([]float64, (n+63)/64)
	parallelFor(n, workers, func(lo, hi int) int64 {
		var u float64
		for i := lo; i < hi; i++ {
			for j := i + 1; j < n; j++ {
				dx := b[j].Pos.X - b[i].Pos.X
				dy := b[j].Pos.Y - b[i].Pos.Y
				u -= g * b[i].Mass * b[j].Mass / math.Sqrt(dx*dx+dy*dy+eps2)
			}
		}
		partial[lo/64] = u
		return 0
	})
	var sum float64
	for _, v := range partial {
		sum += v
	}
	return sum
}

func TotalEnergy(b []dto.Body, g, eps float64, workers int) float64 {
	return Kinetic(b) + Potential(b, g, eps, workers)
}

func Momentum(b []dto.Body) dto.Vec2 {
	var p dto.Vec2
	for _, x := range b {
		p = p.Add(x.Vel.Scale(x.Mass))
	}
	return p
}

func AngularMomentum(b []dto.Body) float64 {
	var l float64
	for _, x := range b {
		l += x.Mass * x.Pos.Cross(x.Vel)
	}
	return l
}
