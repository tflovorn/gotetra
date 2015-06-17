package gotetra

func IterTetras(n, band_index int, Ecache EnergyCache) chan [4]float64 {
	Ets := make(chan [4]float64)
	go iterTetras_worker(n, band_index, Ecache, Ets)
	return Ets
}

func iterTetras_worker(n, band_index int, Ecache EnergyCache, Ets chan [4]float64) {
	// TODO - should vEs be declared before loop for efficency?
	// Placed in loop since it seems like it may hit race condition if
	// outside of inner loop scope.

	// Iterate over tetrahedra.
	subcell_tetras := [6][4]int{[4]int{1, 2, 3, 6}, [4]int{1, 3, 5, 6}, [4]int{3, 5, 6, 7}, [4]int{3, 6, 7, 8}, [4]int{3, 4, 6, 8}, [4]int{2, 3, 4, 6}}
	for k := 0; k < n; k++ {
		for j := 0; j < n; j++ {
			for i := 0; i < n; i++ {
				points := subcell_points(i, j, k)
				for _, sc_t := range subcell_tetras {
					vEs := [4]float64{0.0, 0.0, 0.0, 0.0}
					// Collect band energy at vertices.
					for v, point_index := range sc_t {
						pi, pj, pk := points[point_index-1][0], points[point_index-1][1], points[point_index-1][2]
						vEs[v] = Ecache.EnergyAt(pi, pj, pk, band_index)
					}
					Ets <- vEs
				}
			}
		}
	}
}

func subcell_points(i, j, k int) [8][3]int {
	return [8][3]int{[3]int{i, j, k}, [3]int{i + 1, j, k}, [3]int{i, j + 1, k}, [3]int{i + 1, j + 1, k}, [3]int{i, j, k + 1}, [3]int{i + 1, j, k + 1}, [3]int{i, j + 1, k + 1}, [3]int{i + 1, j + 1, k + 1}}
}

func submesh_index(n, i, j, k int) int {
	return i + j*(n+1) + k*(n+1)*(n+1)
}

func submesh_ijk(n, index int) (int, int, int) {
	i := index % (n + 1)
	j := ((index % ((n + 1) * (n + 1))) - i) / (n + 1)
	k := (index - i - j*(n+1)) / ((n + 1) * (n + 1))
	return i, j, k
}

func submesh_ijk_to_k(n, i, j, k int) [3]float64 {
	step := float64(1 / n)
	return [3]float64{float64(i) * step, float64(j) * step, float64(k) * step}
}

func NumTetra(n int) int {
	// 6 tetrahedra per subcell.
	return 6 * n * n * n
}
