package gotetra

type InputFn func(k [3]float64) []float64

// Returns (total energy, Fermi energy, err).
func SumEnergy(Efn InputFn, n, num_electrons int, R [3][3]float64, use_cache bool) (float64, float64, error) {
	G_order, G_neg := OptimizeGs(R)
	k0 := [3]float64{0.0, 0.0, 0.0}
	num_bands := len(Efn(k0))
	Ecache := NewEnergyCache(n, num_bands, G_order, G_neg, Efn, use_cache)

	E_Fermi, err := FindFermi(n, num_electrons, Ecache)
	if err != nil {
		return 0.0, 0.0, err
	}

	result := 0.0
	c := 0.0
	for i := 0; i < n+1; i++ {
		for j := 0; j < n+1; j++ {
			for k := 0; k < n+1; k++ {
				for band_index := 0; band_index < num_bands; band_index++ {
					this_w := Weight(E_Fermi, n, i, j, k, band_index, Ecache)
					Ekb := Ecache.EnergyAt(i, j, k, band_index)
					contrib := Ekb * this_w
					y := contrib - c
					t := result + y
					c = (t - result) - y
					result = t
				}
			}
		}
	}
	return result, E_Fermi, nil
}
