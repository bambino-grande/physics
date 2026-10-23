package physics

import (
	"math"
	"nbody/physics/dto"
)

const maxDepth = 48

type node struct {
	cx, cy, half float64
	mass         float64
	mx, my       float64
	child        [4]int32
	body         int32
}

func newNode(cx, cy, half float64) node {
	return node{cx: cx, cy: cy, half: half, child: [4]int32{-1, -1, -1, -1}, body: -1}
}

func (n *node) isLeaf() bool {
	return n.child[0] < 0 && n.child[1] < 0 && n.child[2] < 0 && n.child[3] < 0
}

type Tree struct {
	nodes  []node
	bodies []dto.Body
}

func BuildTree(b []dto.Body) *Tree {
	t := &Tree{bodies: b, nodes: make([]node, 0, 2*len(b)+1)}
	if len(b) == 0 {
		t.nodes = append(t.nodes, newNode(0, 0, 1))
		return t
	}
	minX, maxX := b[0].Pos.X, b[0].Pos.X
	minY, maxY := b[0].Pos.Y, b[0].Pos.Y
	for _, p := range b {
		minX, maxX = math.Min(minX, p.Pos.X), math.Max(maxX, p.Pos.X)
		minY, maxY = math.Min(minY, p.Pos.Y), math.Max(maxY, p.Pos.Y)
	}
	half := math.Max(maxX-minX, maxY-minY)/2*1.000001 + 1e-12
	t.nodes = append(t.nodes, newNode((minX+maxX)/2, (minY+maxY)/2, half))
	for i := range b {
		t.insert(int32(i))
	}
	return t
}

func (t *Tree) quadrant(idx int32, x, y float64) int {
	n := &t.nodes[idx]
	q := 0
	if x >= n.cx {
		q |= 1
	}
	if y >= n.cy {
		q |= 2
	}
	return q
}

func (t *Tree) makeChild(parent int32, q int) int32 {
	p := t.nodes[parent]
	h := p.half / 2
	cx, cy := p.cx-h, p.cy-h
	if q&1 != 0 {
		cx = p.cx + h
	}
	if q&2 != 0 {
		cy = p.cy + h
	}
	t.nodes = append(t.nodes, newNode(cx, cy, h))
	id := int32(len(t.nodes) - 1)
	t.nodes[parent].child[q] = id
	return id
}

func (t *Tree) insert(i int32) {
	bi := t.bodies[i]
	idx := int32(0)
	for depth := 0; ; depth++ {
		n := &t.nodes[idx]
		if n.mass == 0 {
			n.body, n.mass, n.mx, n.my = i, bi.Mass, bi.Pos.X, bi.Pos.Y
			return
		}
		total := n.mass + bi.Mass
		n.mx = (n.mx*n.mass + bi.Pos.X*bi.Mass) / total
		n.my = (n.my*n.mass + bi.Pos.Y*bi.Mass) / total
		n.mass = total
		if n.isLeaf() {
			if n.body == -2 || depth >= maxDepth {
				n.body = -2
				return
			}
			j := n.body
			n.body = -1
			pj := t.bodies[j]
			qj := t.quadrant(idx, pj.Pos.X, pj.Pos.Y)
			c := t.makeChild(idx, qj)
			cn := &t.nodes[c]
			cn.body, cn.mass, cn.mx, cn.my = j, pj.Mass, pj.Pos.X, pj.Pos.Y
		}
		q := t.quadrant(idx, bi.Pos.X, bi.Pos.Y)
		if t.nodes[idx].child[q] < 0 {
			t.makeChild(idx, q)
		}
		idx = t.nodes[idx].child[q]
	}
}

func (t *Tree) accel(i int, theta2, g, eps2 float64, stack []int32) (ax, ay float64, inter int64, st []int32) {
	p := t.bodies[i].Pos
	stack = append(stack[:0], 0)
	for len(stack) > 0 {
		idx := stack[len(stack)-1]
		stack = stack[:len(stack)-1]
		n := &t.nodes[idx]
		if n.mass == 0 {
			continue
		}
		leaf := n.isLeaf()
		if leaf && n.body == int32(i) {
			continue
		}
		dx, dy := n.mx-p.X, n.my-p.Y
		if !leaf {
			s := 2 * n.half
			d2 := dx*dx + dy*dy
			inside := math.Abs(p.X-n.cx) <= n.half && math.Abs(p.Y-n.cy) <= n.half
			if inside || s*s >= theta2*d2 {
				for _, c := range n.child {
					if c >= 0 {
						stack = append(stack, c)
					}
				}
				continue
			}
		}
		r2 := dx*dx + dy*dy + eps2
		if r2 == 0 {
			continue
		}
		f := g * n.mass / (r2 * math.Sqrt(r2))
		ax += f * dx
		ay += f * dy
		inter++
	}
	return ax, ay, inter, stack
}

type BarnesHut struct {
	G, Eps  float64
	Theta   float64
	Workers int
}

func (s BarnesHut) Accel(b []dto.Body) dto.Stats {
	t := BuildTree(b)
	eps2, th2 := s.Eps*s.Eps, s.Theta*s.Theta
	inter := parallelFor(len(b), s.Workers, func(lo, hi int) int64 {
		var cnt int64
		stack := make([]int32, 0, 128)
		for i := lo; i < hi; i++ {
			ax, ay, c, st := t.accel(i, th2, s.G, eps2, stack)
			stack = st
			b[i].Acc = dto.Vec2{X: ax, Y: ay}
			cnt += c
		}
		return cnt
	})
	return dto.Stats{Interactions: inter}
}

func (t *Tree) RootMass() float64 { return t.nodes[0].mass }
func (t *Tree) RootCOM() dto.Vec2 { return dto.Vec2{X: t.nodes[0].mx, Y: t.nodes[0].my} }
func (t *Tree) NodeCount() int    { return len(t.nodes) }
