package gotetra

// Return n(E), the total number of states with energy <= E summed over
// all tetrahedra and band indices.
// The calculation of n(E) is implemented as described in BJA94 Appendix A.
//
// TODO doc
// Uses Kahan summation for improved accuracy on dense mesh.
func NumStates(E float64, n int, Ecache EnergyCache, all_bands_at_once bool) float64 {
	num_bands := Ecache.NumBands()
	num_tetra := float64(NumTetra(n))
	result := 0.0
	c := 0.0
	if !all_bands_at_once {
		for band_index := 0; band_index < num_bands; band_index++ {
			for Ets := range IterTetras(n, band_index, Ecache) {
				E1, E2, E3, E4 := Ets[0], Ets[1], Ets[2], Ets[3]
				contrib := NumStatesContrib(E, E1, E2, E3, E4, num_tetra)
				y := contrib - c
				t := result + y
				c = (t - result) - y
				result = t
			}
		}
	} else {
		for Ets := range BandIterTetras(n, Ecache) {
			for band_index := 0; band_index < num_bands; band_index++ {
				E1, E2, E3, E4 := Ets[band_index][0], Ets[band_index][1], Ets[band_index][2], Ets[band_index][3]
				contrib := NumStatesContrib(E, E1, E2, E3, E4, num_tetra)
				y := contrib - c
				t := result + y
				c = (t - result) - y
				result = t
			}
		}
	}
	return result
}

// Return the contribution to the number of states with energy less than
// or equal to E (i.e. the integrated density of states n(E)) from the
// (tetrahedron, band index) pair with the given energies at the vertices.
// The calculation of n(E) is implemented as described in BJA94 Appendix A.
//
// E1, E2, E3, E4 = energies at the vertices of the tetrahedron, in ascending
// order.
//
// num_tetra = total number of tetrahedra in the full Brillouin zone.
// Equal to (volume of tetrahedron) / (volume of full BZ).
func NumStatesContrib(E, E1, E2, E3, E4, num_tetra float64) float64 {
	if E <= E1 {
		return 0.0
	} else if E1 < E && E < E2 {
		return (1.0 / num_tetra) * (E - E1) * (E - E1) * (E - E1) / ((E2 - E1) * (E3 - E1) * (E4 - E1))
	} else if E2 < E && E < E3 {
		fac := (1.0 / num_tetra) / ((E3 - E1) * (E4 - E1))
		esq := (E2-E1)*(E2-E1) + 3*(E2-E1)*(E-E2) + 3*(E-E2)*(E-E2)
		ecub := -(((E3 - E1) + (E4 - E2)) / ((E3 - E2) * (E4 - E2))) * (E - E2) * (E - E2) * (E - E2)
		return fac * (esq + ecub)
	} else if E3 < E && E < E4 {
		return (1.0 / num_tetra) * (1.0 - (E4-E)*(E4-E)*(E4-E)/((E4-E1)*(E4-E2)*(E4-E3)))
	} else {
		// E >= E4
		return (1.0 / num_tetra)
	}
}
