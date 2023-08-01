// Module spectral provides spectral analysis functions for digital signal processing.
module spectral

// segment x segmented into segments of length size with specified noverlap.
// Number of segments returned is (x.len - size) / (size - noverlap) + 1.
pub fn segment(x []f64, size int, noverlap int) [][]f64 {
	stride := size - noverlap
	lx := x.len

	mut segments := 0
	if lx == size {
		segments = 1
	} else if lx > size {
		segments = (x.len - size) / stride + 1
	} else {
		segments = 0
	}

	mut r := [][]f64{len: segments}
	for i, offset := 0, 0; i < segments; i++ {
		r[i] = []f64{len: size}

		for j := 0; j < size; j++ {
			r[i][j] = x[offset + j]
		}

		offset += stride
	}

	return r
}
