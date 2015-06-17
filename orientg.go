package gotetra

// Rearrange the reciprocal lattice vectors such that the Cartesian
// distance from submesh cell point 3 to point 6 (as depicted in BJA94
// Fig. 5) is minimized, in order to minimize interpolation error.
// The reciprocal lattice vectors may be permuted and/or negated.
//
// R = a numpy matrix with rows given by the reciprocal lattice vectors.
//
// Returns two vectors G_order = (o0, o1, o2) and G_neg = (n0, n1, n2).
// The value of G_order is drawn from the permutations of (0, 1, 2) and
// indicates the optimal permutation of the reciprocal lattice vectors,
// when combined with G_neg which indicates whether the corresponding
// reciprocal lattice vector is negated.
// Concretely, R_opt[o0, :] = G_neg[0]*R[0, :] and similarly for R[1, :]
// and R[2, :].
func OptimizeGs(R [3][3]float64) ([3]int, [3]int) {
	permutations := [6][3]int{[3]int{0, 1, 2}, [3]int{0, 2, 1}, [3]int{1, 0, 2}, [3]int{1, 2, 0}, [3]int{2, 0, 1}, [3]int{2, 1, 0}}
	signs := [8][3]int{[3]int{1, 1, 1}, [3]int{1, 1, -1}, [3]int{1, -1, 1}, [3]int{1, -1, -1}, [3]int{-1, 1, 1}, [3]int{-1, 1, -1}, [3]int{-1, -1, 1}, [3]int{-1, -1, -1}}
	G_order := [3]int{0, 0, 0}
	G_neg := [3]int{0, 0, 0}
	opt_36 := 0.0
	init_opt_36 := false
	for _, perm := range permutations {
		for _, sign := range signs {
			this_R := [3][3]float64{[3]float64{0.0, 0.0, 0.0}, [3]float64{0.0, 0.0, 0.0}, [3]float64{0.0, 0.0, 0.0}}
			for i := 0; i < 3; i++ {
				this_R[perm[i]] = vec_mul(R[i], float64(sign[i]))
			}
			mR0 := vec_mul(this_R[0], -1.0)
			R1 := this_R[1]
			mR2 := vec_mul(this_R[2], -1.0)
			k3_to_k6 := vec_norm(vec_add(vec_add(mR0, R1), mR2))
			if !init_opt_36 || k3_to_k6 < opt_36 {
				opt_36 = k3_to_k6
				G_order = perm
				G_neg = sign
				init_opt_36 = true
			}
		}
	}
	return G_order, G_neg
}

// Convert k_opt (a k-point in the represented in the optimal permutation
// of the reciprocal lattice as determined by OptimizeGs) into the
// corresponding k-point in the original reciprocal lattice vectors.
//
// The value of G_order is drawn from the permutations of (0, 1, 2) and
// indicates the optimal permutation of the reciprocal lattice vectors,
// when combined with G_neg which indicates whether the corresponding
// reciprocal lattice vector is negated.
// Concretely, R_opt[o0, :] = G_neg[0]*R[0, :] and similarly for R[1, :]
// and R[2, :].
func Get_k_Orig(k_opt [3]float64, G_order, G_neg [3]int) [3]float64 {
	k_orig := [3]float64{0.0, 0.0, 0.0}
	for i := 0; i < 3; i++ {
		k_orig[i] = k_opt[G_order[i]] * float64(G_neg[i])
	}
	return k_orig
}

// Return the value of Ropt corresponding to G_order and G_neg.
//
// The value of G_order is drawn from the permutations of (0, 1, 2) and
// indicates the optimal permutation of the reciprocal lattice vectors,
// when combined with G_neg which indicates whether the corresponding
// reciprocal lattice vector is negated.
// Concretely, R_opt[o0, :] = G_neg[0]*R[0, :] and similarly for R[1, :]
// and R[2, :].
func GetRopt(R [3][3]float64, G_order, G_neg [3]int) [3][3]float64 {
	R_opt := [3][3]float64{[3]float64{0.0, 0.0, 0.0}, [3]float64{0.0, 0.0, 0.0}, [3]float64{0.0, 0.0, 0.0}}
	for i := 0; i < 3; i++ {
		R_opt[G_order[i]] = vec_mul(R[i], float64(G_neg[i]))
	}
	return R_opt
}
