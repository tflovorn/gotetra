package gotetra

// Return n(E), the total number of states with energy <= E summed over
// all tetrahedra and band indices.
// The calculation of n(E) is implemented as described in BJA94 Appendix A.
//
// TODO doc
func NumStates(E float64, n int, Efn InputFn, Ecache *EnergyCache) float64 {
	result := 0.0
	num_bands := Ecache.NumBands()
	num_tetra := NumTetra(n)
	for band_index := 0; band_index < num_bands; band_index++ {
		for Ets := range IterTetras(n, band_index, Efn, Ecache) {
			// TODO - make sure that returning Ets from IterTetras
			// isn't allocating. Don't expect it to be - since
			// Ets is fixed length, should go on the stack
			// (created as literal [4]float64{E1, E2, E3, E4}).
			E1, E2, E3, E4 := Ets[0], Ets[1], Ets[2], Ets[3]
			result += NumStatesContrib(E, Ets, num_tetra)
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
func NumStatesContrib(E, E1, E2, E3, E4, num_tetra float64) {
	if E <= E1 {
		return 0.0
	} else if E1 < E && E < E2 {
		return (1 / num_tetra) * (E - E1) * (E - E1) * (E - E1) / ((E2 - E1) * (E3 - E1) * (E4 - E1))
	} else if E2 < E && E < E3 {
		fac = (1 / num_tetra) / ((E3 - E1) * (E4 - E1))
		esq = (E2-E1)*(E2-E1) + 3*(E2-E1)*(E-E2) + 3*(E-E2)*(E-E2)
		ecub = -(((E3 - E1) + (E4 - E2)) / ((E3 - E2) * (E4 - E2))) * (E - E2) * (E - E2) * (E - E2)
		return fac * (esq + ecub)
	} else if E3 < E && E < E4 {
		return (1 / num_tetra) * (1 - (E4-E)*(E4-E)*(E4-E)/((E4-E1)*(E4-E2)*(E4-E3)))
	} else {
		// E >= E4
		return (1 / num_tetra)
	}
}
