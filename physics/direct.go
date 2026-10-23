package physics

import (
	"math"
	"nbody/physics/dto"
)

type Direct struct {
	G, Eps  float64
	Workers int
}

func (d Direct) Accel(b []dto.Body) dto.Stats {
	eps2 := d.Eps * d.Eps
	n := len(b)
	inter := parallelFor(n, d.Workers, func(lo, hi int) int64 {
		var cnt int64
		for i := lo; i < hi; i++ {
			var ax, ay float64
			pi := b[i].Pos
			for j := 0; j < n; j++ {
				if j == i {
					continue
				}
				dx := b[j].Pos.X - pi.X
				dy := b[j].Pos.Y - pi.Y
				r2 := dx*dx + dy*dy + eps2
				if r2 == 0 {
					continue
				}
				f := d.G * b[j].Mass / (r2 * math.Sqrt(r2))
				ax += f * dx
				ay += f * dy
				cnt++
			}
			b[i].Acc = dto.Vec2{X: ax, Y: ay}
		}
		return cnt
	})
	return dto.Stats{Interactions: inter}
}
