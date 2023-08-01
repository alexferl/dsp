module fft

import benchmark
import math
import math.complex
import utils

const sqrt2_2 = math.sqrt2 / 2

struct FFTTest {
	inp []f64
	out []complex.Complex
}

const fft_tests = [
	// vfmt off
	// impulse responses
	FFTTest{
		[
			f64(1)
		],
		[
			complex.complex(1, 0)
		]
	},
	FFTTest{
		[
			f64(1), 0
		],
		[
			complex.complex(1, 0),
			complex.complex(1, 0)
		]
	},
	FFTTest{
		[
			f64(1), 0, 0, 0
		],
		[
			complex.complex(1, 0),
			complex.complex(1, 0),
			complex.complex(1, 0),
			complex.complex(1, 0)
		]
	},
	FFTTest{
		[
			f64(1), 0, 0, 0, 0, 0, 0, 0
		],
		[
			complex.complex(1, 0),
			complex.complex(1, 0),
			complex.complex(1, 0),
			complex.complex(1, 0),
			complex.complex(1, 0),
			complex.complex(1, 0),
			complex.complex(1, 0),
			complex.complex(1, 0)
		]
	},
	// shifted impulse response
	FFTTest{
		[
			f64(0), 1
		],
		[
			complex.complex(1, 0)
			complex.complex(-1, 0)
		]
	},
	FFTTest{
		[
			f64(0), 1, 0, 0
		],
		[
			complex.complex(1, 0)
			complex.complex(0, -1)
			complex.complex(-1, 0)
			complex.complex(0, 1)
		]
	},
	FFTTest{
		[
			f64(0), 1, 0, 0, 0, 0, 0, 0
		],
		[
			complex.complex(1, 0)
			complex.complex(sqrt2_2, -sqrt2_2)
			complex.complex(0, -1)
			complex.complex(-sqrt2_2, -sqrt2_2)
			complex.complex(-1, 0)
			complex.complex(-sqrt2_2, sqrt2_2)
			complex.complex(0, 1)
			complex.complex(sqrt2_2, sqrt2_2)
		]
	},
	// other
	FFTTest{
		[
			f64(1), 2, 3, 4
		],
		[
			complex.complex(10, 0)
			complex.complex(-2, 2)
			complex.complex(-2, 0)
			complex.complex(-2, -2)
		]
	},
	FFTTest{
		[
			f64(1), 3, 5, 7
		],
		[
			complex.complex(16, 0)
			complex.complex(-4, 4)
			complex.complex(-4, 0)
			complex.complex(-4, -4)
		]
	},
	FFTTest{
		[
			f64(1), 2, 3, 4, 5, 6, 7, 8
		],
		[
			complex.complex(36, 0),
			complex.complex(-4, 9.65685425),
			complex.complex(-4, 4),
			complex.complex(-4, 1.65685425),
			complex.complex(-4, 0),
			complex.complex(-4, -1.65685425),
			complex.complex(-4, -4),
			complex.complex(-4, -9.65685425)
		]
	},
	// non power of 2 lengths
	FFTTest{
		[
			f64(1), 0, 0, 0, 0
		],
		[
			complex.complex(1, 0)
			complex.complex(1, 0)
			complex.complex(1, 0)
			complex.complex(1, 0)
			complex.complex(1, 0)
		]
	},
	FFTTest{
		[
			f64(1), 2, 3
		],
		[
			complex.complex(6, 0)
			complex.complex(-1.5, 0.8660254)
			complex.complex(-1.5, -0.8660254)
		]
	},
	FFTTest{
		[
			f64(1), 1, 1
		],
		[
			complex.complex(3, 0)
			complex.complex(0, 0)
			complex.complex(0, 0)
		]
	},
	// vfmt on
]

fn test_fft() {
	tol := 1e-8
	for t in fft.fft_tests {
		v := fft_real(t.inp)
		assert utils.equal_approx_complex(v, t.out, tol), 'input: ${t.inp}\n output: ${v}\n expected: ${t.out}'

		vi := ifft(t.out)
		assert utils.equal_approx_complex(vi, utils.to_complex(t.inp), tol), 'input: ${t.out}\n output: ${vi}\n expected: ${utils.to_complex(t.inp)}'
	}
}

struct FFT2Test {
	inp [][]f64
	out [][]complex.Complex
}

const fft2_tests = [
	// vfmt off
	FFT2Test{
		[
			[f64(1), 2, 3]
			[f64(3), 4, 5]
		],
		[
			[complex.complex(18, 0), complex.complex(-3, 1.73205081), complex.complex(-3, -1.73205081)]
			[complex.complex(-6, 0), complex.complex(0, 0), complex.complex(0, 0)]
		]
	},
	FFT2Test{
		[
			[f64(0.1), 0.2, 0.3, 0.4, 0.5]
			[f64(1), 2, 3, 4, 5]
			[f64(3), 2, 1, 0, -1]
		],
		[
			[complex.complex(21.5, 0), complex.complex(-0.25, 0.34409548), complex.complex(-0.25, 0.08122992), complex.complex(-0.25, -0.08122992), complex.complex(-0.25, -0.34409548)]
			[complex.complex(-8.5, -8.66025404), complex.complex(5.70990854, 4.6742225), complex.complex(1.15694356, 4.41135694), complex.complex(-1.65694356, 4.24889709), complex.complex(-6.20990854, 3.98603154)]
			[complex.complex(-8.5, 8.66025404), complex.complex(-6.20990854, -3.98603154), complex.complex(-1.65694356, -4.24889709), complex.complex(1.15694356, -4.41135694), complex.complex(5.70990854, -4.6742225)]
		]
	}
	// vfmt on
]

fn test_fft2() {
	tol := 1e-8
	for t in fft.fft2_tests {
		v := fft2_real(t.inp)
		assert utils.equal_approx_2d_complex(v, t.out, tol), 'fft2 error\ninput: ${t.inp}\noutput: ${v}\nexpected: ${t.out}'

		vi := ifft2_real(t.out)
		assert utils.equal_approx_2d_complex(vi, utils.to_complex_2d(t.inp), tol), 'ifft2 error\ninput: ${t.out}\noutput: ${vi}\nexpected: ${utils.to_complex_2d(t.inp)}'
	}
}

struct FFTNTest {
	inp []f64
	dim []int
	out []complex.Complex
}

const fftn_tests = [
	// vfmt off
	FFTNTest{
			[
				f64(4), 2, 3, 8, 5, 6, 7, 2, 13, 24, 13, 17
			]
			[
				int(2), 2, 3
			]
			[
				complex.complex(104, 0), complex.complex(12.5, 14.72243186), complex.complex(12.5, -14.72243186)
				complex.complex(-42, 0), complex.complex(-10.5, 6.06217783), complex.complex(-10.5, -6.06217783)
				complex.complex(-48, 0), complex.complex(-4.5, -11.25833025), complex.complex(-4.5, 11.25833025)
				complex.complex(22, 0), complex.complex(8.5, -6.06217783), complex.complex(8.5, 6.06217783)
			]
	}
	// vfmt on
]

fn test_fftn() {
	tol := 1e-8
	for t in fft.fftn_tests {
		m := utils.new_matrix(utils.to_complex(t.inp), t.dim)
		o := utils.new_matrix(t.out, t.dim)
		v := fftn(m)
		assert v.equal_approx(o, tol), 'fftn error\ninput: ${m}\noutput: ${v}\nexpected: ${o}'

		vi := ifftn(o)
		assert vi.equal_approx(m, tol), 'ifftn error\ninput: ${o}\noutput: ${vi}\nexpected: ${m}'
	}
}

struct ReverseBitsTest {
	inp u32
	sz  u32
	out u32
}

const reverse_bits_tests = [
	ReverseBitsTest{0, 1, 0},
	ReverseBitsTest{1, 2, 2},
	ReverseBitsTest{1, 4, 8},
	ReverseBitsTest{2, 4, 4},
	ReverseBitsTest{3, 4, 12},
]

fn test_reverse_bits() {
	for t in fft.reverse_bits_tests {
		v := reverse_bits(t.inp, t.sz)
		assert v == t.out, 'input: ${t.inp} size: ${t.sz} output: ${v} expected: ${t.out}'
	}
}

fn test_benchmark_fft() {
	mut bmark := benchmark.new_benchmark()
	bmark.set_total_expected_steps(10)
	for i := 0; i < bmark.nexpected_steps; i++ {
		big_n := 1 << 20
		mut a := []complex.Complex{len: big_n}
		for j := 0; j < big_n; j++ {
			a[j] = complex.complex(f64(j) / f64(big_n), 0)
		}

		mut r := new_radix2()
		r.get_factors(big_n)

		bmark.step()
		fft(a)
		bmark.ok()
	}
	bmark.stop()
	bmark.measure(@FN)
}
