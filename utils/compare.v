module utils

import math
import math.complex

// equal_approx returns true when the slices have equal lengths and
// all element pairs have an absolute tolerance less than tol.
pub fn equal_approx(a []f64, b []f64, tol f64) bool {
	if a.len != b.len {
		return false
	}

	for i, c in a {
		if !equal_within_abs(c, b[i], tol) {
			return false
		}
	}

	return true
}

// equal_approx_complex returns true when the slices have equal lengths and
// all element pairs have an absolute tolerance less than tol.
pub fn equal_approx_complex(a []complex.Complex, b []complex.Complex, tol f64) bool {
	if a.len != b.len {
		return false
	}

	for i, c in a {
		if !equal_within_abs_complex(c, b[i], tol) {
			return false
		}
	}

	return true
}

// equal_approx_2d returns true when the 2d slices have equal lengths and
// all element pairs have an absolute tolerance less than tol.
pub fn equal_approx_2d(a [][]f64, b [][]f64, tol f64) bool {
	if a.len != b.len {
		return false
	}

	for i, c in a {
		if !equal_approx(c, b[i], tol) {
			return false
		}
	}

	return true
}

// equal_approx_2d_complex returns true when the 2d slices have equal lengths and
// all element pairs have an absolute tolerance less than tol.
pub fn equal_approx_2d_complex(a [][]complex.Complex, b [][]complex.Complex, tol f64) bool {
	if a.len != b.len {
		return false
	}

	for i, c in a {
		if !equal_approx_complex(c, b[i], tol) {
			return false
		}
	}

	return true
}

// equal_within_abs returns true when a and b have an absolute difference
// not greater than tol.
pub fn equal_within_abs(a f64, b f64, tol f64) bool {
	return math.abs(a - b) <= tol || math.abs(1 - a / b) <= tol
}

// equal_within_abs_complex returns true when a and b have an absolute difference
// not greater than tol.
pub fn equal_within_abs_complex(a complex.Complex, b complex.Complex, tol f64) bool {
	r_a := a.re
	r_b := b.re
	i_a := a.im
	i_b := b.im
	return equal_within_abs(r_a, r_b, tol) && equal_within_abs(i_a, i_b, tol)
}
