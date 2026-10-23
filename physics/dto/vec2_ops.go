package dto

import "math"

func (a Vec2) Add(b Vec2) Vec2      { return Vec2{a.X + b.X, a.Y + b.Y} }
func (a Vec2) Sub(b Vec2) Vec2      { return Vec2{a.X - b.X, a.Y - b.Y} }
func (a Vec2) Scale(s float64) Vec2 { return Vec2{a.X * s, a.Y * s} }
func (a Vec2) Dot(b Vec2) float64   { return a.X*b.X + a.Y*b.Y }
func (a Vec2) Norm() float64        { return math.Sqrt(a.X*a.X + a.Y*a.Y) }
func (a Vec2) Cross(b Vec2) float64 { return a.X*b.Y - a.Y*b.X }
