package gotetra

import (
	"math"
)

func vec_mul(v [3]float64, s float64) [3]float64 {
	u := [3]float64{v[0] * s, v[1] * s, v[2] * s}
	return u
}

func vec_add(v, u [3]float64) [3]float64 {
	w := [3]float64{v[0] + u[0], v[1] + u[1], v[2] + u[2]}
	return w
}

func vec_norm(v [3]float64) float64 {
	norm := math.Sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
	return norm
}
