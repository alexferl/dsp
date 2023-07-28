module fft

import math
import math.complex
import runtime
import sync

[noinit]
struct Radix2 {
mut:
	mu      sync.RwMutex
	factors map[int][]complex.Complex = {
		4: [complex.complex(1, 0), complex.complex(0, 1), complex.complex(-1, 0),
			complex.complex(0, 1)]
	}
}

fn new_radix2() Radix2 {
	return Radix2{}
}

fn (mut r Radix2) get_factors(input_len int) []complex.Complex {
	r.mu.@rlock()

	if r.has_factors(input_len) {
		defer {
			r.mu.runlock()
		}
		return r.factors[input_len]
	}

	r.mu.runlock()
	r.mu.@lock()
	defer {
		r.mu.unlock()
	}

	if !r.has_factors(input_len) {
		for i, p := i32(8), 4; i <= input_len; i, p = i << 1, i {
			_ := r.factors[i] or {
				r.factors[i] = []complex.Complex{len: int(i)}

				for n, j := 0, 0; n < i; n, j = n + 2, j + 1 {
					r.factors[i][n] = r.factors[p][j]
				}

				for n := 1; n < i; n += 2 {
					sin, cos := math.sincos(-2 * math.pi / f64(i) * f64(n))
					r.factors[i][n] = complex.complex(cos, sin)
				}

				return r.factors[input_len]
			}
		}
	}

	return r.factors[input_len]
}

fn (r Radix2) has_factors(idx int) bool {
	_ := r.factors[idx] or { return false }
	return true
}

struct FFTWork {
	start  int
	end    int
}

const worker_pool_size = 0

fn (mut r Radix2) fft(x []complex.Complex) []complex.Complex {
	lx := x.len
	factors := r.get_factors(lx)

	mut t := []complex.Complex{len: lx}
	mut rd := reorder_data(x)

	mut blocks, mut stage, mut s_2 := int(0), int(0), int(0)

	mut ch := sync.new_channel[FFTWork](u32(lx))
	mut wg := sync.new_waitgroup()

	mut num_workers := fft.worker_pool_size
	if num_workers == 0 {
		num_workers = runtime.nr_jobs()
	}

	mut idx_diff := lx / num_workers
	if idx_diff < 2 {
		idx_diff = 2
	}

	worker := fn [mut wg, mut ch, stage, s_2, rd, factors, blocks, mut t]() {
		for {
			mut work := FFTWork{}
			if !ch.pop(&work) {
				break
			}

			for nb := work.start; nb < work.end; nb += stage {
				if stage != 2 {
					for j := 0; j < s_2; j++ {
						idx := j + nb
						idx2 := idx + s_2
						ridx := rd[idx]
						w_n := rd[idx2] * factors[blocks * j]
						t[idx] = ridx + w_n
						t[idx2] = ridx - w_n
					}
				} else {
					n1 := nb + 1
					rn := rd[nb]
					rn1 := rd[n1]
					t[nb] = rn + rn1
					t[n1] = rn - rn1
				}
			}
		}
		wg.done()
	}

	for i := 0; i < num_workers; i++ {
		spawn worker()
	}
	defer { ch.close() }

	for stage = 2; stage <= lx; stage <<= 1 {
		blocks = lx / stage
		s_2 = stage / 2

		for start, end := 0, stage; true; {
			if end - start >= idx_diff || end == lx {
				wg.add(1)
				work := FFTWork{start, end}
				ch.push(&work)

				if end == lx {
					break
				}

				start = end
			}

			end += stage
		}
		// ch.close()
		wg.wait()
		rd = t.clone()
		t = rd.clone()
	}

	return rd
}

fn reorder_data(x []complex.Complex) []complex.Complex {
	lx := u32(x.len)
	mut r := []complex.Complex{len: int(lx)}
	s := log2(lx)

	mut n := u32(0)
	for ; n < lx; n++ {
		r[reverse_bits(n, u32(s))] = x[n]
	}

	return r
}

// log2 returns the log base 2 of v
// from: http://graphics.stanford.edu/~seander/bithacks.html#IntegerLogObvious
fn log2(v_ u32) u32 {
	mut v := v_
	mut r := u32(0)

	// TODO: change to v >>= 1 when it works
	for v = v >> 1; v != 0; v >>= 1 {
		r++
	}

	return r
}

// reverse_bits returns the first s bits of v in reverse order
// from: http://graphics.stanford.edu/~seander/bithacks.html#BitReverseObvious
fn reverse_bits(v_ u32, s_ u32) u32 {
	mut v := v_
	mut s := s_
	mut r := u32(0)

	// Since we aren't reversing all the bits in v (just the first s bits),
	// we only need the first bit of v instead of a full copy.
	r = v & 1
	s--

	// TODO: change to v >>= 1 when it works
	for v = v >> 1; v != 0; v >>= 1 {
		r <<= 1
		r |= v & 1
		s--
	}

	return r << s
}
