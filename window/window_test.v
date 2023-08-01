module window

import utils

struct WindowTest {
	inp      int
	hamming  []f64
	hann     []f64
	bartlett []f64
	flat_top []f64
	blackman []f64
}

const window_tests = [
	// vfmt off
	WindowTest{
		1,
		[f64(1)],
		[f64(1)],
		[f64(1)],
		[f64(1)],
		[f64(1)]
	},
	WindowTest{
		5,
		[f64(0.08), 0.54, 1, 0.54, 0.08],
		[f64(0), 0.5, 1, 0.5, 0],
		[f64(0), 0.5, 1,0.5, 0],
		[f64(-0.0004210510000000013), -0.05473684000000003, 1, -0.05473684000000003, -0.0004210510000000013],
		[f64(0), 0.34, 1, 0.34, 0]
	},
	WindowTest{
		10,
		[f64(0.08), 0.18761956, 0.46012184, 0.77, 0.97225861, 0.97225861, 0.77, 0.46012184, 0.18761956, 0.08],
		[f64(0), 0.116977778440511, 0.413175911166535, 0.75, 0.969846310392954, 0.969846310392954, 0.75, 0.413175911166535, 0.116977778440511, 0],
		[f64(0), 0.222222222222222, 0.444444444444444, 0.666666666666667, 0.888888888888889, 0.888888888888889, 0.666666666666667, 0.444444444444444, 0.222222222222222, 0],
		[f64(-0.000421051000000), -0.020172031509486, -0.070199042063189, 0.198210530000000, 0.862476344072674, 0.862476344072674, 0.198210530000000, -0.070199042063189, -0.020172031509486, -0.000421051000000],
		[f64(0), 0.0508696327, 0.258000502, 0.63, 0.951129866, 0.951129866, 0.63, 0.258000502, 0.0508696327, 0]},
	// vfmt on
]

fn test_window_functions() {
	tol := 1e-8
	for t in window.window_tests {
		mut o := hamming(t.inp)
		assert utils.equal_approx(o, t.hamming, tol), 'hamming error\ninput: ${t.inp}\noutput: ${o}\nexpected: ${t.hamming}'

		o = hann(t.inp)
		assert utils.equal_approx(o, t.hann, tol), 'hann error\ninput: ${t.inp}\noutput: ${o}\nexpected: ${t.hann}'

		o = bartlett(t.inp)
		assert utils.equal_approx(o, t.bartlett, tol), 'bartlett error\ninput: ${t.inp}\noutput: ${o}\nexpected: ${t.bartlett}'

		o = rectangular(t.inp)
		apply(mut o, hamming)
		assert utils.equal_approx(o, t.hamming, tol), 'apply error\noutput: ${o}\nexpected: ${t.hamming}'

		o = flat_top(t.inp)
		assert utils.equal_approx(o, t.flat_top, tol), 'flat_top error\ninput: ${t.inp}\noutput: ${o}\nexpected: ${t.flat_top}'

		o = blackman(t.inp)
		assert utils.equal_approx(o, t.blackman, tol), 'blackman error\ninput: ${t.inp}\noutput: ${o}\nexpected: ${t.flat_top}'
	}
}
