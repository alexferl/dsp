module spectral

import utils

struct SegmentTest {
	size     int
	noverlap int
	out      [][]f64
}

const segment_tests = [
	SegmentTest{4, 0, [[f64(1), 2, 3, 4], [f64(5), 6, 7, 8]]},
	SegmentTest{4, 1, [[f64(1), 2, 3, 4], [f64(5), 6, 7, 8], [f64(7), 8, 9, 10]]},
	SegmentTest{4, 2, [[f64(1), 2, 3, 4], [f64(3), 4, 5, 6], [f64(5), 6, 7, 8],
		[f64(7), 8, 9, 10]]},
]

fn test_segment() {
	tol := 1e-8
	x := [f64(1), 2, 3, 4, 5, 6, 7, 8, 9, 10]

	for t in spectral.segment_tests {
		o := segment(x, t.size, t.noverlap)
		utils.equal_approx_2d(o, t.out, tol), 'segment error: output: ${o} expected: ${t.out}'
	}
}
