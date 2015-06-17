package gotetra

import (
	"testing"
)

func TestSum_SimpleInsulator(t *testing.T) {
	n := 4
	a := 1.0
	R := cubicR(a)
	num_bands := 2
	num_electrons := 1
	t0 := 1.0
	E0 := 6.0
	deltaE := 14.0 // need deltaE > bandwidth = 12t for insulator
	Efn := func(k [3]float64) []float64 {
		return simplebands_Efn(k, num_bands, t0, E0, deltaE)
	}

	Esum, E_Fermi, err := SumEnergy(Efn, n, num_electrons, R)
	if err != nil {
		t.Fatal(err)
	}
	println("Sum", Esum)
	println("E_Fermi", E_Fermi)
}
