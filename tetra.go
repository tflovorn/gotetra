package gotetra

func NumTetra(n int) int {
	// 6 tetrahedra per subcell.
	return 6 * n * n * n
}

func IterTetras(n, band_index int, Ecache EnergyCache) chan [4]float64 {
	Ets := make(chan [4]float64)
	go iterTetras_worker(n, band_index, Ecache, Ets, nil)
	return Ets
}

func iterTetras_worker(n, band_index int, Ecache EnergyCache, Ets_chan chan [4]float64, ks_chan chan [4][3]int) {
	// TODO - should vEs be declared before loop for efficency?
	// Placed in loop since it seems like it may hit race condition if
	// outside of inner loop scope.
	tetras := subcell_tetras()
	// Iterate over tetrahedra.
	for k := 0; k < n; k++ {
		for j := 0; j < n; j++ {
			for i := 0; i < n; i++ {
				points := subcell_points(i, j, k)
				for _, sc_t := range tetras {
					vEs := [4]float64{0.0, 0.0, 0.0, 0.0}
					ks := [4][3]int{[3]int{0, 0, 0}, [3]int{0, 0, 0}, [3]int{0, 0, 0}, [3]int{0, 0, 0}}
					// Collect band energy at vertices.
					for v, point_index := range sc_t {
						pi, pj, pk := points[point_index-1][0], points[point_index-1][1], points[point_index-1][2]
						ks[v] = [3]int{pi, pj, pk}
						vEs[v] = Ecache.EnergyAt(pi, pj, pk, band_index)
					}
					// Sort and yield vertex energies.
					sortEsKs(&vEs, &ks)
					Ets_chan <- vEs
					// Yield associated ks if requested.
					if ks_chan != nil {
						ks_chan <- ks
					}
				}
			}
		}
	}
	close(Ets_chan)
	if ks_chan != nil {
		close(ks_chan)
	}
}

func IterTetrasKs(n, band_index int, Ecache EnergyCache) (chan [4]float64, chan [4][3]int) {
	Ets_chan, ks_chan := make(chan [4]float64), make(chan [4][3]int)
	go iterTetras_worker(n, band_index, Ecache, Ets_chan, ks_chan)
	return Ets_chan, ks_chan
}

func sortEsKs(Es *[4]float64, ks *[4][3]int) {
	// insertion sort
	for i := 1; i < 4; i++ {
		j := i
		for j > 0 && Es[j-1] > Es[j] {
			// swap j, j-1
			Es[j], Es[j-1] = Es[j-1], Es[j]
			ks[j], ks[j-1] = ks[j-1], ks[j]
			j -= 1
		}
	}
}

func IterTetrasAround(n, i, j, k, band_index int, Ecache EnergyCache) (chan [4]float64, chan [4][3]int) {
	Ets_chan, ks_chan := make(chan [4]float64), make(chan [4][3]int)
	go iterTetrasAround_worker(n, i, j, k, band_index, Ecache, Ets_chan, ks_chan)
	return Ets_chan, ks_chan
}

func iterTetrasAround_worker(n, i, j, k, band_index int, Ecache EnergyCache, Ets_chan chan [4]float64, ks_chan chan [4][3]int) {
	tetras := subcell_tetras()
	subcells := subcells_around_ijk(n, i, j, k)
	for _, sc := range subcells {
		sc_i, sc_j, sc_k := sc[0], sc[1], sc[2]
		sc_points := subcell_points(sc_i, sc_j, sc_k)
		for _, sc_t := range tetras {
			vEs := [4]float64{0.0, 0.0, 0.0, 0.0}
			ks := [4][3]int{[3]int{0, 0, 0}, [3]int{0, 0, 0}, [3]int{0, 0, 0}, [3]int{0, 0, 0}}
			// Collect band energy at vertices.
			for v, point_index := range sc_t {
				pi, pj, pk := sc_points[point_index-1][0], sc_points[point_index-1][1], sc_points[point_index-1][2]
				ks[v] = [3]int{pi, pj, pk}
				vEs[v] = Ecache.EnergyAt(pi, pj, pk, band_index)
			}
			// Sort and yield vertex energies.
			sortEsKs(&vEs, &ks)
			Ets_chan <- vEs
			// Yield associated ks if requested.
			if ks_chan != nil {
				ks_chan <- ks
			}
		}
	}
	close(Ets_chan)
	if ks_chan != nil {
		close(ks_chan)
	}
}

func subcell_tetras() [6][4]int {
	return [6][4]int{[4]int{1, 2, 3, 6}, [4]int{1, 3, 5, 6}, [4]int{3, 5, 6, 7}, [4]int{3, 6, 7, 8}, [4]int{3, 4, 6, 8}, [4]int{2, 3, 4, 6}}
}

func subcell_points(i, j, k int) [8][3]int {
	return [8][3]int{[3]int{i, j, k}, [3]int{i + 1, j, k}, [3]int{i, j + 1, k}, [3]int{i + 1, j + 1, k}, [3]int{i, j, k + 1}, [3]int{i + 1, j, k + 1}, [3]int{i, j + 1, k + 1}, [3]int{i + 1, j + 1, k + 1}}
}

func submesh_index(n, i, j, k int) int {
	return i + j*(n+1) + k*(n+1)*(n+1)
}

func submesh_ijk(n, index int) (int, int, int) {
	i := index % (n + 1)
	j := round(float64((index%((n+1)*(n+1)))-i) / float64(n+1))
	k := round(float64(index-i-j*(n+1)) / float64((n+1)*(n+1)))
	return i, j, k
}

func submesh_ijk_to_k(n, i, j, k int) [3]float64 {
	step := 1.0 / float64(n)
	return [3]float64{float64(i) * step, float64(j) * step, float64(k) * step}
}

func subcells_around_ijk(n, i, j, k int) [][3]int {
	subcells := [][3]int{}
	// For each subcell, point (i, j, k) has the number in BJA94 Fig. 5
	// corresponding to the subcell number.
	// Subcell #1
	if i != n && j != n && k != n {
		subcells = append(subcells, [3]int{i, j, k})
	}
	// Subcell #2
	if i != 0 && j != n && k != n {
		subcells = append(subcells, [3]int{i - 1, j, k})
	}
	// Subcell #3
	if i != n && j != 0 && k != n {
		subcells = append(subcells, [3]int{i, j - 1, k})
	}
	// Subcell #4
	if i != 0 && j != 0 && k != n {
		subcells = append(subcells, [3]int{i - 1, j - 1, k})
	}
	// Subcell #5
	if i != n && j != n && k != 0 {
		subcells = append(subcells, [3]int{i, j, k - 1})
	}
	// Subcell #6
	if i != 0 && j != n && k != 0 {
		subcells = append(subcells, [3]int{i - 1, j, k - 1})
	}
	// Subcell #7
	if i != n && j != 0 && k != 0 {
		subcells = append(subcells, [3]int{i, j - 1, k - 1})
	}
	// Subcell #8
	if i != 0 && j != 0 && k != 0 {
		subcells = append(subcells, [3]int{i - 1, j - 1, k - 1})
	}
	return subcells
}
