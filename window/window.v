// Module window provides window functions for digital signal processing.
module window

import math

// apply applies the window f to x.
pub fn apply(mut x []f64, f fn (int) []f64) {
	for i, w in f(x.len) {
		x[i] *= w
	}
}

// rectangular returns an L-point rectangular window (all values are 1).
pub fn rectangular(big_l int) []f64 {
	mut r := []f64{len: big_l}

	for i, _ in r {
		r[i] = 1
	}

	return r
}

// hamming returns an L-point symmetric Hamming window.
// Reference: http://www.mathworks.com/help/signal/ref/hamming.html
pub fn hamming(big_l int) []f64 {
	mut r := []f64{len: big_l}

	if big_l == 1 {
		r[0] = 1
	} else {
		big_n := big_l - 1
		coef := math.pi * 2 / f64(big_n)
		for n := 0; n <= big_n; n++ {
			r[n] = 0.54 - 0.46 * math.cos(coef * f64(n))
		}
	}

	return r
}

// hann returns an L-point Hann window.
// Reference: http://www.mathworks.com/help/signal/ref/hann.html
pub fn hann(big_l int) []f64 {
	mut r := []f64{len: big_l}

	if big_l == 1 {
		r[0] = 1
	} else {
		big_n := big_l - 1
		coef := 2 * math.pi / f64(big_n)
		for n := 0; n <= big_n; n++ {
			r[n] = 0.5 * (1 - math.cos(coef * f64(n)))
		}
	}

	return r
}

// bartlett returns an L-point Bartlett window.
// Reference: http://www.mathworks.com/help/signal/ref/bartlett.html
pub fn bartlett(big_l int) []f64 {
	mut r := []f64{len: big_l}

	if big_l == 1 {
		r[0] = 1
	} else {
		big_n := big_l - 1
		coef := 2 / f64(big_n)
		mut n := 0
		for ; n <= big_n / 2; n++ {
			r[n] = coef * f64(n)
		}
		for ; n <= big_n; n++ {
			r[n] = 2 - coef * f64(n)
		}
	}

	return r
}

// FlatTop returns an L-point flat top window.
// Reference: http://www.mathworks.com/help/signal/ref/flattopwin.html
pub fn flat_top(big_l int) []f64 {
	alpha0 := f64(0.21557895)
	alpha1 := f64(0.41663158)
	alpha2 := f64(0.277263158)
	alpha3 := f64(0.083578947)
	alpha4 := f64(0.006947368)

	mut r := []f64{len: big_l}

	if big_l == 1 {
		r[0] = 1
		return r
	}

	big_n := big_l - 1
	coef := 2 * math.pi / f64(big_n)

	for n := 0; n <= big_n; n++ {
		factor := f64(n) * coef

		term0 := alpha0
		term1 := alpha1 * math.cos(factor)
		term2 := alpha2 * math.cos(2 * factor)
		term3 := alpha3 * math.cos(3 * factor)
		term4 := alpha4 * math.cos(4 * factor)

		r[n] = term0 - term1 + term2 - term3 + term4
	}

	return r
}

// blackman returns an L-point Blackman window
// Reference: http://www.mathworks.com/help/signal/ref/blackman.html
pub fn blackman(big_l int) []f64 {
	mut r := []f64{len: big_l}

	if big_l == 1 {
		r[0] = 1
	} else {
		big_n := big_l - 1
		for n := 0; n <= big_n; n++ {
			term0 := 0.42
			term1 := -0.5 * math.cos(2 * math.pi * f64(n) / f64(big_n))
			term2 := 0.08 * math.cos(4 * math.pi * f64(n) / f64(big_n))
			r[n] = term0 + term1 + term2
		}
	}

	return r
}
