module utils

import math.complex

struct SegmentTest {
	segs     int
	noverlap f64
	slices   [][]int
}

const segment_tests = [
	// vfmt off
	SegmentTest{3, .5, [[int(0), 8], [int(4), 12], [int(8), 16]]},
	// vfmt on
]

fn test_segment() {
	tol := 1e-8
	mut x := []complex.Complex{len: 16}
	for n, _ in x {
		x[n] = complex.complex(f64(n), 0)
	}

	for _, st in utils.segment_tests {
		v := segment(x, st.segs, st.noverlap)
		mut s := [][]complex.Complex{len: st.segs}
		for i, sl in st.slices {
			s[i] = x[sl[0]..sl[1]]
		}

		assert equal_approx_2d_complex(v, s, tol), 'segment error: expected: ${s} output: ${v}'
	}
}
