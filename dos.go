package gotetra

// Return the contribution to the density of states at energy E from the
// (tetrahedron, band index) pair with the given energies at the vertices.
// The calculation of D_T(E) is implemented as described in BJA94 Appendix C.
//
// E1, E2, E3, E4 = energies at the vertices of the tetrahedron, in ascending
// order.
//
// num_tetra = total number of tetrahedra in the full Brillouin zone.
// Equal to (volume of tetrahedron) / (volume of full BZ).
func DosContrib(E, E1, E2, E3, E4, num_tetra float64) float64 {
	if E <= E1 {
		return 0.0
	} else if E1 < E && E < E2 {
		return (1 / num_tetra) * 3 * (E - E1) * (E - E1) / ((E2 - E1) * (E3 - E1) * (E4 - E1))
	} else if E2 < E && E < E3 {
		fac := (1 / num_tetra) / ((E3 - E1) * (E4 - E1))
		elin := 3*(E2-E1) + 6*(E-E2)
		esq := -3 * (((E3 - E1) + (E4 - E2)) / ((E3 - E2) * (E4 - E2))) * (E - E2) * (E - E2)
		return fac * (elin + esq)
	} else if E3 < E && E < E4 {
		return (1 / num_tetra) * 3 * (E4 - E) * (E4 - E) / ((E4 - E1) * (E4 - E2) * (E4 - E3))
	} else {
		// E >= E4
		return 0.0
	}
}
