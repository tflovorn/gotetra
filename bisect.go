package gotetra

import (
	"errors"
	"math"
)

// Root-finding by bisection.
// Returns a value x such that a <= x <= b and |errfn(x)| < toly.
// Following the implementation in Quantum Espresso, no tolerance for the
// root bracket interval size is used (unlike Scipy's bisect).
// TODO - is this better? Why? Expect that it may be due to insulators -
// large range of x's where errfn(x) = 0 (in gap).
func Bisect(errfn func(float64) float64, a, b, toly float64, maxiter int) (float64, error) {
	val_a := errfn(a)
	val_b := errfn(b)
	if val_a == 0.0 {
		return a, nil
	} else if val_b == 0.0 {
		return b, nil
	}
	if (val_a < 0.0 && val_b < 0.0) || (val_a > 0.0 && val_b > 0.0) {
		return 0.0, errors.New("a and b do not bracket root in Bisect()")
	}

	low, high := a, b // val_a < 0
	if val_b < 0.0 {
		low, high = b, a
	}
	for iter := 0; iter < maxiter; iter++ {
		mid := (low + high) / 2
		val_mid := errfn(mid)
		if math.Abs(val_mid) < toly {
			return mid, nil
		} else if val_mid > 0.0 {
			// go toward low
			high = mid
		} else {
			// go toward high
			low = mid
		}
	}
	return 0.0, errors.New("maximum iterations reached in Bisect()")
}
