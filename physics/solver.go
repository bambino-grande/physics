package physics

import (
	"nbody/physics/dto"
	"runtime"
	"sync/atomic"
)

type Solver interface {
	Accel(b []dto.Body) dto.Stats
}

func NumCPU() int { return runtime.NumCPU() }

func parallelFor(n, workers int, fn func(lo, hi int) int64) int64 {
	if workers <= 1 || n < 128 {
		return fn(0, n)
	}
	const chunk = 64
	var next, total int64
	done := make(chan struct{}, workers)
	for w := 0; w < workers; w++ {
		go func() {
			var local int64
			for {
				lo := int(atomic.AddInt64(&next, chunk)) - chunk
				if lo >= n {
					break
				}
				hi := lo + chunk
				if hi > n {
					hi = n
				}
				local += fn(lo, hi)
			}
			atomic.AddInt64(&total, local)
			done <- struct{}{}
		}()
	}
	for w := 0; w < workers; w++ {
		<-done
	}
	return total
}
