package gotetra

type EnergyCache interface {
	EnergyAt(i, j, k, band_index int) float64
	NumBands() int
}

type energyCacheRam struct {
	n         int
	num_bands int
	G_order   [3]int
	G_neg     [3]int
	Efn       InputFn
	Eks       [][]float64
}

func NewEnergyCache(n, num_bands int, G_order, G_neg [3]int, Efn InputFn) EnergyCache {
	// Default to RAM cache for now.
	// TODO - check n to choose RAM cache or disk cache.
	ec := new(energyCacheRam)
	nks := (n + 1) * (n + 1) * (n + 1)
	ec.n = n
	ec.num_bands = num_bands
	ec.G_order = G_order
	ec.G_neg = G_neg
	ec.Efn = Efn
	// Only need to pre-allocate Eks for k-points.
	// Will set (k, band) eigenvalue list after k-point is queried.
	Eks := make([][]float64, nks)
	ec.Eks = Eks
	return ec
}

func (ec *energyCacheRam) EnergyAt(i, j, k, band_index int) float64 {
	k_opt := submesh_ijk_to_k(ec.n, i, j, k)
	k_index := submesh_index(ec.n, i, j, k)
	if ec.Eks[k_index] == nil {
		// Already queried this (i, k, k).
		return ec.Eks[k_index][band_index]
	}
	// Haven't seen this (i, j, k) before; set it.
	k_orig := Get_k_Orig(k_opt, ec.G_order, ec.G_neg)
	Es := ec.Efn(k_orig)
	ec.Eks[k_index] = Es
	return Es[band_index]
}

func (ec *energyCacheRam) NumBands() int {
	return ec.num_bands
}
