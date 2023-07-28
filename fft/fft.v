module fft

import arrays
import math.complex
import utils

// fft_real returns the forward FFT of the real-valued slice.
pub fn fft_real(x []f64) []complex.Complex {
	return fft(utils.to_complex(x))
}

// ifft_real returns the inverse FFT of the real-valued slice.
pub fn ifft_real(x []f64) []complex.Complex {
	return ifft(utils.to_complex(x))
}

pub fn ifft(x []complex.Complex) []complex.Complex {
	lx := x.len
	mut r := []complex.Complex{len: lx}

	// Reverse inputs, which is calculated with modulo N, hence x[0] as an outlier
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

pub fn fft(x []complex.Complex) []complex.Complex {
	lx := x.len

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

// fft2real returns the 2-dimensional, forward FFT of the real-valued matrix.
pub fn fft2real(x [][]f64) [][]complex.Complex {
	return fft2(utils.to_complex_2(x))
}

// fft2 returns the 2-dimensional, forward FFT of the complex-valued matrix.
pub fn fft2(x [][]complex.Complex) [][]complex.Complex {
	return compute_fft2(x, fft)
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
