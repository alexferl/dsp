module utils

import arrays
import math
import math.complex

// to_complex returns the complex equivalent of the real-valued slice.
pub fn to_complex(x []f64) []complex.Complex {
	mut y := []complex.Complex{len: x.len}
	for n, v in x {
		y[n] = complex.complex(v, 0)
	}
	return y
}

// is_power_of_2 returns true if x is a power of 2, else false.
pub fn is_power_of_2(x int) bool {
	return x & (x - 1) == 0
}

// next_power_of_2 returns the next power of 2 >= x.
pub fn next_power_of_2(x int) int {
	if is_power_of_2(x) {
		return x
	}

	return int(math.pow(2, math.ceil(math.log2(f64(x)))))
}

// zero_pad returns x with zeros appended to the end to the specified length.
// If len(x) >= length, x is returned, otherwise a new array is returned.
pub fn zero_pad[T](x []T, length int) []T {
	if x.len >= length {
		return x
	}

	mut r := []T{len: length}
	arrays.copy[T](mut r, x)
	return r
}

// zero_pad_2 returns zero_pad of x, with the length as the next power of 2 >= len(x).
pub fn zero_pad_2(x []complex.Complex) []complex.Complex {
	return zero_pad[complex.Complex](x, next_power_of_2(x.len))
}

// to_complex_2 returns the complex equivalent of the real-valued matrix.
pub fn to_complex_2(x [][]f64) [][]complex.Complex {
	mut y := [][]complex.Complex{len: x.len}
	for n, v in x {
		y[n] = to_complex(v)
	}
	return y
}

// segment returns segs equal-length slices that are segments of x with noverlap% of overlap.
// The returned slices are not copies of x, but slices into it.
// Trailing entries in x that cannot be included in the equal-length segments are discarded.
// noverlap is a percentage, thus 0 <= noverlap <= 1, and noverlap = 0.5 is 50% overlap.
pub fn segment(x []complex.Complex, segs int, noverlap f64) [][]complex.Complex {
	lx := x.len

	mut overlap := 0
	mut length := 0
	mut step := 0
	mut tot := 0
	for length = lx; length > 0; length-- {
		overlap = int(f64(length) * noverlap)
		tot = segs * (length - overlap) + overlap
		if tot <= lx {
			step = length - overlap
			break
		}
	}

	if length == 0 {
		panic('too many segments')
	}

	mut r := [][]complex.Complex{len: segs}
	mut s := 0
	for n, _ in r {
		r[n] = x[s..s + length]
		s += step
	}

	return r
}
