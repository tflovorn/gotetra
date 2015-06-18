package gotetra

func Weight(E_Fermi float64, n, i, j, k, band_index int, Ecache EnergyCache) float64 {
	result := 0.0
	c := 0.0
	num_tetra := float64(NumTetra(n))
	Ets_chan, ks_chan := IterTetrasAround(n, i, j, k, band_index, Ecache)
	// TODO - will calculate weight contributions more times than
	// necessary. Could solve with weight contributions cache.
	// (check if (ks, num_bands) weight has been evaluated already).

	// Iterate over all tetrahedra. For each one, if it includes a
	// contribution to the weight at this k-point, add that contribution.
	// Ignore the other tetrahedra, and ignore contributions to other
	// k-points sharing tetrahedra with this one.
	for Ets := range Ets_chan {
		ks := <-ks_chan
		for index_k, this_kval := range ks {
			this_i, this_j, this_k := this_kval[0], this_kval[1], this_kval[2]
			if this_i == i && this_j == j && this_k == k {
				E1, E2, E3, E4 := Ets[0], Ets[1], Ets[2], Ets[3]
				contribs := WeightContrib(E_Fermi, E1, E2, E3, E4, num_tetra)
				my_contrib := contribs[index_k]
				y := my_contrib - c
				t := result + y
				c = (t - result) - y
				result = t
			}
		}
	}
	return result
}

func BandWeights(E_Fermi float64, n, i, j, k int, Ecache EnergyCache) ([]float64, []float64) {
	num_bands := Ecache.NumBands()
	result := make([]float64, num_bands)
	c := make([]float64, num_bands)
	num_tetra := float64(NumTetra(n))
	Ets_chan, ks_chan := BandIterTetrasAround(n, i, j, k, Ecache)
	this_Es_set := false
	this_Es := make([]float64, num_bands)
	// TODO - will calculate weight contributions more times than
	// necessary. Could solve with weight contributions cache.
	// (check if (ks, num_bands) weight has been evaluated already).

	// Iterate over all tetrahedra. For each one, if it includes a
	// contribution to the weight at this k-point, add that contribution.
	// Ignore the other tetrahedra, and ignore contributions to other
	// k-points sharing tetrahedra with this one.
	for Ets := range Ets_chan {
		ks := <-ks_chan
		for band_index := 0; band_index < num_bands; band_index++ {
			for index_k, this_kval := range ks[band_index] {
				this_i, this_j, this_k := this_kval[0], this_kval[1], this_kval[2]
				if this_i == i && this_j == j && this_k == k {
					E1, E2, E3, E4 := Ets[band_index][0], Ets[band_index][1], Ets[band_index][2], Ets[band_index][3]
					contribs := WeightContrib(E_Fermi, E1, E2, E3, E4, num_tetra)
					my_contrib := contribs[index_k]
					y := my_contrib - c[band_index]
					t := result[band_index] + y
					c[band_index] = (t - result[band_index]) - y
					result[band_index] = t
					if !this_Es_set {
						for bi := 0; bi < num_bands; bi++ {
							this_Es[bi] = Ets[bi][index_k]
						}
						this_Es_set = true
					}
				}
			}
		}
	}
	return result, this_Es
}

// Return the specified tetrahedron's contribution to the integration
// weights at the k-points of the tetrahedron's vertices; i.e. return
// a list with elements w_{bandIndex, kN, tetra}, where the elements of the
// returned list range over kN values in the order specified by tetra.
// The calculation of the tetrahedron contribution to w_{nj} is implemented
// as described in BJA94 Appendix B and Section V.
//
// E_Fermi = Fermi energy of the system.
//
// TODO docs
func WeightContrib(E_Fermi, E1, E2, E3, E4, num_tetra float64) [4]float64 {
	ws := [4]float64{0.0, 0.0, 0.0, 0.0}
	qnt := 1.0 / (4.0 * num_tetra)
	if E_Fermi <= E1 {
		for i := 0; i < 4; i++ {
			// redundant here - leaving in for clarity
			ws[i] = 0.0
		}
	} else if E1 < E_Fermi && E_Fermi < E2 {
		C := qnt * (E_Fermi - E1) * (E_Fermi - E1) * (E_Fermi - E1) / ((E2 - E1) * (E3 - E1) * (E4 - E1))
		ws[0] = C * (4.0 - (E_Fermi-E1)*(1.0/(E2-E1)+1.0/(E3-E1)+1.0/(E4-E1)))
		ws[1] = C * (E_Fermi - E1) / (E2 - E1)
		ws[2] = C * (E_Fermi - E1) / (E3 - E1)
		ws[3] = C * (E_Fermi - E1) / (E4 - E1)
	} else if E2 < E_Fermi && E_Fermi < E3 {
		C1, C2, C3 := cs_23(E_Fermi, E1, E2, E3, E4, qnt)
		ws[0] = C1 + (C1+C2)*(E3-E_Fermi)/(E3-E1) + (C1+C2+C3)*(E4-E_Fermi)/(E4-E1)
		ws[1] = C1 + C2 + C3 + (C2+C3)*(E3-E_Fermi)/(E3-E2) + C3*(E4-E_Fermi)/(E4-E2)
		ws[2] = (C1+C2)*(E_Fermi-E1)/(E3-E1) + (C2+C3)*(E_Fermi-E2)/(E3-E2)
		ws[3] = (C1+C2+C3)*(E_Fermi-E1)/(E4-E1) + C3*(E_Fermi-E2)/(E4-E2)
	} else if E3 < E_Fermi && E_Fermi < E4 {
		C := qnt * (E4 - E_Fermi) * (E4 - E_Fermi) * (E4 - E_Fermi) / ((E4 - E1) * (E4 - E2) * (E4 - E3))
		ws[0] = qnt - C*(E4-E_Fermi)/(E4-E1)
		ws[1] = qnt - C*(E4-E_Fermi)/(E4-E2)
		ws[2] = qnt - C*(E4-E_Fermi)/(E4-E3)
		ws[3] = qnt - C*(4.0-(1.0/(E4-E1)+1.0/(E4-E2)+1.0/(E4-E3))*(E4-E_Fermi))
	} else {
		// E_Fermi >= E4
		for i := 0; i < 4; i++ {
			ws[i] = qnt
		}
	}

	dws := curvatureCorrection(E_Fermi, E1, E2, E3, E4, num_tetra)
	for i := 0; i < 4; i++ {
		ws[i] += dws[i]
	}

	return ws
}

func cs_23(E_Fermi, E1, E2, E3, E4, qnt float64) (float64, float64, float64) {
	C1 := qnt * (E_Fermi - E1) * (E_Fermi - E1) / ((E4 - E1) * (E3 - E1))

	C2_num := qnt * (E_Fermi - E1) * (E_Fermi - E2) * (E3 - E_Fermi)
	C2_denom := (E4 - E1) * (E3 - E2) * (E3 - E1)
	C2 := C2_num / C2_denom

	C3_num := qnt * (E_Fermi - E2) * (E_Fermi - E2) * (E4 - E_Fermi)
	C3_denom := (E4 - E2) * (E3 - E2) * (E4 - E1)
	C3 := C3_num / C3_denom

	return C1, C2, C3
}

// Return a list of the curvature corrections to the k-point weight
// contributions from the specified tetrahedron. The band energies at the
// vertices of the tetrahedron are given in sorted order by Es; the returned
// corrections are given in the same order.
// TODO fix doc
func curvatureCorrection(E_Fermi, E1, E2, E3, E4, num_tetra float64) [4]float64 {
	D_T := DosContrib(E_Fermi, E1, E2, E3, E4, num_tetra)
	dws := [4]float64{0.0, 0.0, 0.0, 0.0}
	Es := [4]float64{E1, E2, E3, E4}
	sumEs := E1 + E2 + E3 + E4
	for i := 0; i < 4; i++ {
		dws[i] = (D_T / 40.0) * (sumEs - 4.0*Es[i])
	}
	return dws
}
