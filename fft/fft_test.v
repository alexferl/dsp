module fft

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
	FFT2Test{},
]

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
