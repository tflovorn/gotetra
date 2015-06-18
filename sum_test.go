package gotetra

import (
	"math"
	"testing"
)

func TestSum_SimpleInsulator(t *testing.T) {
	all_bands_at_once := true
	use_cache := true
	n := 32
	a := 1.0
	R := cubicR(a)
	num_bands := 10
	num_electrons := 1
	t0 := 1.0
	E0 := 6.0
	deltaE := 14.0 // need deltaE > bandwidth = 12t for insulator
	Efn := func(k [3]float64) []float64 {
		return simplebands_Efn(k, num_bands, t0, E0, deltaE)
	}

	Esum, E_Fermi, err := SumEnergy(Efn, n, num_electrons, R, all_bands_at_once, use_cache)
	if err != nil {
		t.Fatal(err)
	}

	expected_sum := E0
	eps := 1e-9
	if math.Abs(Esum-expected_sum) > eps {
		t.Fatalf("Incorrect energy sum")
	}
	if E_Fermi < E0+6.0*t0 || E_Fermi > E0+deltaE-6.0*t0 {
		t.Fatalf("Incorrect Fermi energy position")
	}
}
