package gotetra

import (
	"sort"
)

type EnergyCache interface {
	EnergyAt(i, j, k, band_index int) float64
	NumBands() int
	MinE() float64
	MaxE() float64
}

type energyCacheRam struct {
	use_cache bool
	n         int
	num_bands int
	G_order   [3]int
	G_neg     [3]int
	Efn       InputFn
	Eks       [][]float64
}

func NewEnergyCache(n, num_bands int, G_order, G_neg [3]int, Efn InputFn, use_cache bool) EnergyCache {
	// Default to RAM cache for now.
	// TODO - check n to choose RAM cache or disk cache.
	ec := new(energyCacheRam)
	nks := (n + 1) * (n + 1) * (n + 1)
	ec.n = n
	ec.num_bands = num_bands
	ec.G_order = G_order
	ec.G_neg = G_neg
	ec.Efn = Efn
	ec.use_cache = use_cache
	if use_cache {
		// Only need to pre-allocate Eks for k-points.
		// Will set (k, band) eigenvalue list after k-point is queried.
		Eks := make([][]float64, nks)
		ec.Eks = Eks
	} else {
		Eks := [][]float64{}
		ec.Eks = Eks
	}
	return ec
}

func (ec *energyCacheRam) EnergyAt(i, j, k, band_index int) float64 {
	k_opt := submesh_ijk_to_k(ec.n, i, j, k)
	k_index := submesh_index(ec.n, i, j, k)
	if ec.use_cache && ec.Eks[k_index] != nil {
		// Already queried this (i, k, k).
		return ec.Eks[k_index][band_index]
	}
	// Haven't seen this (i, j, k) before; set it.
	k_orig := Get_k_Orig(k_opt, ec.G_order, ec.G_neg)
	Es := ec.Efn(k_orig)
	if !sort.Float64sAreSorted(Es) {
		panic("Got unsorted values from Efn in EnergyAt()")
	}
	if ec.use_cache {
		ec.Eks[k_index] = Es
	}
	return Es[band_index]
}

func (ec *energyCacheRam) NumBands() int {
	return ec.num_bands
}

func (ec *energyCacheRam) MinE() float64 {
	minval := 0.0
	minval_init := false
	for i := 0; i < ec.n+1; i++ {
		for j := 0; j < ec.n+1; j++ {
			for k := 0; k < ec.n+1; k++ {
				for b := 0; b < ec.num_bands; b++ {
					E := ec.EnergyAt(i, j, k, b)
					if !minval_init || E < minval {
						minval = E
						minval_init = true
					}
				}
			}
		}
	}
	return minval
}

func (ec *energyCacheRam) MaxE() float64 {
	maxval := 0.0
	maxval_init := false
	for i := 0; i < ec.n+1; i++ {
		for j := 0; j < ec.n+1; j++ {
			for k := 0; k < ec.n+1; k++ {
				for b := 0; b < ec.num_bands; b++ {
					E := ec.EnergyAt(i, j, k, b)
					if !maxval_init || E > maxval {
						maxval = E
						maxval_init = true
					}
				}
			}
		}
	}
	return maxval
}
