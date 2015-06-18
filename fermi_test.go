package gotetra

import (
	"math"
	"testing"
)

// TODO - will need a matrix inverse for general D.
// gonum? Inverse(m Matrix) (*Dense, error); (m *Dense) Row(i int) --> []float64
func cubicR(a float64) [3][3]float64 {
	Ra := [3]float64{2.0 * math.Pi / a, 0.0, 0.0}
	Rb := [3]float64{0.0, 2.0 * math.Pi / a, 0.0}
	Rc := [3]float64{0.0, 0.0, 2.0 * math.Pi / a}
	R := [3][3]float64{Ra, Rb, Rc}
	return R
}

func simplebands_Efn(k [3]float64, num_bands int, t, E0, deltaE float64) []float64 {
	Eks := make([]float64, num_bands)
	tk := -2 * t * (math.Cos(2.0*math.Pi*k[0]) + math.Cos(2.0*math.Pi*k[1]) + math.Cos(2.0*math.Pi*k[2]))
	for i := 0; i < num_bands; i++ {
		E0_band := E0 + float64(i)*deltaE
		Eks[i] = E0_band + tk
	}
	return Eks
}

func TestFermi_SimpleInsulator(t *testing.T) {
	all_bands_at_once := true
	use_cache := true
	n := 8
	a := 1.0
	R := cubicR(a)
	G_order, G_neg := OptimizeGs(R)
	num_bands := 2
	num_electrons := 1.0
	t0 := 1.0
	E0 := 6.0
	deltaE := 14.0 // need deltaE > bandwidth = 12t for insulator
	Efn := func(k [3]float64) []float64 {
		return simplebands_Efn(k, num_bands, t0, E0, deltaE)
	}
	Ecache := NewEnergyCache(n, num_bands, G_order, G_neg, Efn, use_cache)

	fermi, err := FindFermi(n, num_electrons, Ecache, all_bands_at_once)
	if err != nil {
		t.Fatal(err)
	}
	if fermi < E0+6.0*t0 || fermi > E0+deltaE-6.0*t0 {
		t.Fatalf("Incorrect Fermi energy position")
	}
}
