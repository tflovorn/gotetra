package gotetra

type InputFn func(k [3]float64) []float64

// Returns (total energy, Fermi energy, err).
func SumEnergy(Efn InputFn, n int, num_electrons float64, R [3][3]float64, all_bands_at_once, use_cache bool) (float64, float64, error) {
	G_order, G_neg := OptimizeGs(R)
	k0 := [3]float64{0.0, 0.0, 0.0}
	num_bands := len(Efn(k0))
	Ecache := NewEnergyCache(n, num_bands, G_order, G_neg, Efn, use_cache)

	E_Fermi, err := FindFermi(n, num_electrons, Ecache, all_bands_at_once)
	if err != nil {
		return 0.0, 0.0, err
	}

	result := 0.0
	if !all_bands_at_once {
		result = individualBandSum(n, num_bands, E_Fermi, Ecache)
	} else {
		result = simultaneousBandSum(n, num_bands, E_Fermi, Ecache)
	}
	return result, E_Fermi, nil
}

func individualBandSum(n, num_bands int, E_Fermi float64, Ecache EnergyCache) float64 {
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
	return result
}

func simultaneousBandSum(n, num_bands int, E_Fermi float64, Ecache EnergyCache) float64 {
	result := 0.0
	c := 0.0
	for i := 0; i < n+1; i++ {
		for j := 0; j < n+1; j++ {
			for k := 0; k < n+1; k++ {
				this_ws, this_Es := BandWeights(E_Fermi, n, i, j, k, Ecache)
				for band_index := 0; band_index < num_bands; band_index++ {
					contrib := this_Es[band_index] * this_ws[band_index]
					y := contrib - c
					t := result + y
					c = (t - result) - y
					result = t
				}
			}
		}
	}
	return result
}
