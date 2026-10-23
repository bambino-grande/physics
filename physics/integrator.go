package physics

import "nbody/physics/dto"

func Step(b []dto.Body, dt float64, s Solver) dto.Stats {
	h := dt / 2
	for i := range b {
		b[i].Vel.X += b[i].Acc.X * h
		b[i].Vel.Y += b[i].Acc.Y * h
		b[i].Pos.X += b[i].Vel.X * dt
		b[i].Pos.Y += b[i].Vel.Y * dt
	}
	st := s.Accel(b)
	for i := range b {
		b[i].Vel.X += b[i].Acc.X * h
		b[i].Vel.Y += b[i].Acc.Y * h
	}
	return st
}

func Init(b []dto.Body, s Solver) dto.Stats { return s.Accel(b) }
