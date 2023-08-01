module fft

import math
import math.complex
import sync
import utils

[noinit]
struct Bluestein {
mut:
	mu          sync.RwMutex
	factors     map[int][]complex.Complex
	inv_factors map[int][]complex.Complex
}

fn new_bluestein() Bluestein {
	return Bluestein{}
}

fn (mut b Bluestein) get_factors(input_len int) ([]complex.Complex, []complex.Complex) {
	b.mu.@rlock()

	if b.has_factors(input_len) {
		defer {
			b.mu.runlock()
		}
		return b.factors[input_len], b.inv_factors[input_len]
	}

	b.mu.runlock()
	b.mu.@lock()
	defer {
		b.mu.unlock()
	}

	if !b.has_factors(input_len) {
		b.factors[input_len] = []complex.Complex{len: input_len}
		b.inv_factors[input_len] = []complex.Complex{len: input_len}

		mut sin := f64(0)
		mut cos := f64(0)
		for i := 0; i < input_len; i++ {
			if i == 0 {
				sin, cos = 0, 1
			} else {
				sin, cos = math.sincos(math.pi / f64(input_len) * f64(i * i))
			}
			b.factors[input_len][i] = complex.complex(cos, sin)
			b.inv_factors[input_len][i] = complex.complex(cos, -sin)
		}
	}

	return b.factors[input_len], b.inv_factors[input_len]
}

fn (b Bluestein) has_factors(idx int) bool {
	_ := b.factors[idx] or { return false }
	return true
}

fn (mut b Bluestein) fft(x []complex.Complex) []complex.Complex {
	lx := x.len
	mut a := utils.zero_pad[complex.Complex](x, utils.next_power_of_two(lx * 2 - 1))
	la := a.len
	factors, inv_factors := b.get_factors(lx)

	for n, v in x {
		a[n] = v * inv_factors[n]
	}

	mut b1 := []complex.Complex{len: la}
	for i := 0; i < lx; i++ {
		b1[i] = factors[i]

		if i != 0 {
			b1[la - i] = factors[i]
		}
	}

	mut r := convolve(a, b1)

	for i := 0; i < lx; i++ {
		r[i] *= inv_factors[i]
	}

	return r[..lx]
}
