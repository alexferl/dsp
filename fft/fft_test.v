module fft

import math.complex
import math
import utils

const sqrt2_2 = math.sqrt2 / 2

struct FFTTest {
	inp []f64
	out []complex.Complex
}

const fft_tests = [
	FFTTest{[f64(1)], [complex.complex(1, 0)]},
	FFTTest{[f64(1), 0], [complex.complex(1, 0), complex.complex(1, 0)]},
	// FFTTest{[f64(1), 0, 0, 0], [complex.complex(1, 0), complex.complex(1, 0), complex.complex(1, 0), complex.complex(1, 0)]},
]

fn test_fft() {
	tol := 1e-8
	for t in fft.fft_tests {
		v := fft_real(t.inp)
		assert utils.equal_approx_complex(v, t.out, tol), 'input: ${t.inp} output: ${v} expected: ${t.out}'

		vi := ifft(t.out)
		assert utils.equal_approx_complex(vi, utils.to_complex(t.inp), tol), 'input: ${t.out} output: ${vi} expected: ${utils.to_complex(t.inp)}'
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
