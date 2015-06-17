package gotetra

// Returns the Fermi energy E_F, at which the integrated number of
// states n(E_F) = num_electrons.
// Assumes that Ecache has access to all the electronic states of the system;
// e.g. for a collinear spin-polarized system, Ecache must not contain only
// up-spins or only down-spins but instead contain all states of both types.
//
// TODO doc
func FindFermi(n, num_electrons int, Ecache EnergyCache) (float64, error) {
	statecount_error := func(E float64) float64 {
		count := NumStates(E, n, Ecache)
		return count - float64(num_electrons)
	}
	emin := Ecache.MinE()
	emax := Ecache.MaxE()
	toln := 1e-10
	maxiter := 300
	E_Fermi, err := Bisect(statecount_error, emin, emax, toln, maxiter)
	if err != nil {
		return 0.0, err
	}
	return E_Fermi, nil
}
