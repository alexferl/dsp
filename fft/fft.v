module fft

import arrays
import math.complex
import utils

pub fn fft(x []complex.Complex) []complex.Complex {
	lx := x.len

	// TODO: non-hack handling length <= 1 cases
	if lx <= 1 {
		mut r := []complex.Complex{len: lx}
		arrays.copy[complex.Complex](mut r, x)
		return r
	}

	mut r := new_radix2()
	if utils.is_power_of_2(lx) {
		return r.fft(x)
	}

	mut b := new_bluestein()
	return b.fft(x)
}

// fft_real returns the forward FFT of the real-valued slice.
pub fn fft_real(x []f64) []complex.Complex {
	return fft(utils.to_complex(x))
}

pub fn ifft(x []complex.Complex) []complex.Complex {
	lx := x.len
	mut r := []complex.Complex{len: lx}

	// Reverse inputs, which is calculated with modulo n, hence x[0] as an outlier
	r[0] = x[0]
	for i := 1; i < lx; i++ {
		r[i] = x[lx - i]
	}

	r = fft(r)

	big_n := complex.complex(f64(lx), 0)
	for n, _ in r {
		r[n] /= big_n
	}

	return r
}

// ifft_real returns the inverse FFT of the real-valued slice.
pub fn ifft_real(x []f64) []complex.Complex {
	return ifft(utils.to_complex(x))
}

// convolve returns the convolution of x ∗ y.
pub fn convolve(x []complex.Complex, y []complex.Complex) []complex.Complex {
	if x.len != y.len {
		panic('arrays not of equal size"')
	}

	fft_x := fft(x)
	fft_y := fft(y)

	mut r := []complex.Complex{len: x.len}
	for i := 0; i < r.len; i++ {
		r[i] = fft_x[i] * fft_y[i]
	}

	return ifft(r)
}

// fft2 returns the 2-dimensional, forward FFT of the complex-valued matrix.
pub fn fft2(x [][]complex.Complex) [][]complex.Complex {
	return compute_fft2(x, fft)
}

// fft2real returns the 2-dimensional, forward FFT of the real-valued matrix.
pub fn fft2real(x [][]f64) [][]complex.Complex {
	return fft2(utils.to_complex_2(x))
}

// ifft2real returns the 2-dimensional, inverse FFT of the real-valued matrix.
pub fn ifft2real(x [][]complex.Complex) [][]complex.Complex {
	return compute_fft2(x, ifft)
}

fn compute_fft2(x [][]complex.Complex, f fn ([]complex.Complex) []complex.Complex) [][]complex.Complex {
	rows := x.len
	if rows == 0 {
		panic('empty input array')
	}

	cols := x[0].len
	mut r := [][]complex.Complex{len: rows}
	for i := 0; i < rows; i++ {
		if x[i].len != cols {
			panic('ragged input array')
		}
		r[i] = []complex.Complex{len: cols}
	}

	for i := 0; i < cols; i++ {
		mut t := []complex.Complex{len: rows}
		for j := 0; j < rows; j++ {
			t[j] = x[j][i]
		}

		for n, v in f(t) {
			r[n][i] = v
		}
	}

	for n, v in r {
		r[n] = f(v)
	}

	return r
}

// fftn returns the forward FFT of the matrix m, computed in all N dimensions.
pub fn fftn(m utils.Matrix) utils.Matrix {
	return compute_fttn(m, fft)
}

// ifftn returns the forward FFT of the matrix m, computed in all N dimensions.
pub fn ifftn(m utils.Matrix) utils.Matrix {
	return compute_fttn(m, ifft)
}

fn compute_fttn(m utils.Matrix, f fn ([]complex.Complex) []complex.Complex) utils.Matrix {
	mut dims := m.dimensions()
	mut t := m.copy()
	mut r := utils.new_empty_matrix(dims)

	for n in dims {
		dims[n] -= 1
	}

	for n in dims {
		mut d := []int{len: dims.len}
		arrays.copy[int](mut d, dims)
		d[n] = -1

		for {
			r.set_dim(f(t.dim(d)), d)

			if !decr_dim(mut d, dims) {
				break
			}
		}

		r, t = t, r
	}

	return t
}

// decr_dim decrements an element of x by 1, skipping all -1s, and wrapping up to d.
// If a value is 0, it will be reset to its corresponding value in d, and will carry one from the next non -1 value to the right.
// Returns true if decremented, else false.
fn decr_dim(mut x []int, d []int) bool {
	for n, v in x {
		if v == -1 {
			continue
		} else if v == 0 {
			mut i := n
			// find the next element to decrement
			for ; i < x.len; i++ {
				if x[i] == -1 {
					continue
				} else if x[i] == 0 {
					x[i] = d[i]
				} else {
					x[i] -= 1
					return true
				}
			}

			// no decrement
			return false
		} else {
			x[n] -= 1
			return true
		}
	}

	return false
}
