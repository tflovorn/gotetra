package gotetra

import (
	"math"
	"testing"
)

func TestBisectCubic(t *testing.T) {
	cubic := func(x float64) float64 {
		return x*x*x + 1
	}
	toly := 1e-12
	maxiter := 300
	a := -5.0
	b := 5.0
	root, err := Bisect(cubic, a, b, toly, maxiter)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(cubic(root)) > toly {
		t.Fatalf("Bisect() did not reach tolerance")
	}
	expected := -1.0
	if math.Abs(root-expected) > 1e-6 {
		t.Fatalf("Bisect found incorrect root")
	}
}
