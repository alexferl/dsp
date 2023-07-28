module utils

import math.complex
import arrays

// Matrix is a multidimensional matrix of arbitrary size and dimension.
// It cannot be resized after creation. Arrays in any axis can be set or fetched.
pub struct Matrix {
mut:
	list    []complex.Complex
	dims    []int
	offsets []int
}

// new_matrix returns a new Matrix populated with x having dimensions dims.
// For example, to create a 3-dimensional Matrix with 2 components, 3 rows, and 4 columns:
//   new_matrix([]complex.Complex {
//     1, 2, 3, 4,
//     5, 6, 7, 8,
//     9, 0, 1, 2,
//
//     3, 4, 5, 6,
//     7, 8, 9, 0,
//     4, 3, 2, 1},
//   []int {2, 3, 4})
pub fn new_matrix(x []complex.Complex, dims []int) Matrix {
	mut length := 1
	mut offsets := []int{len: dims.len}

	for i := dims.len - 1; i >= 0; i-- {
		if dims[i] < 1 {
			panic('invalid dimensions')
		}

		offsets[i] = length
		length *= dims[i]
	}

	if x.len != length {
		panic('incorrect dimensions')
	}

	mut dc := []int{len: dims.len}
	arrays.copy[int](mut dc, dims)
	return Matrix{
		list: x
		dims: dc
		offsets: offsets
	}
}

// new_matrix_2d is a helper function to convert a 2d array to a matrix.
pub fn new_matrix_2d(x [][]complex.Complex) Matrix {
	dims := [x.len, x[0].len]
	mut r := []complex.Complex{len: dims[0] * dims[1]}
	for n, v in x {
		if v.len != dims[1] {
			panic('ragged array')
		}

		arrays.copy[complex.Complex](mut r[n * dims[1]..(n + 1) * dims[1]], v)
	}

	return new_matrix(r, dims)
}

// new_empty_matrix creates an empty Matrix with given dimensions.
pub fn new_empty_matrix(dims []int) Matrix {
	mut x := 1
	for v in dims {
		x *= v
	}

	return new_matrix([]complex.Complex{len: x}, dims)
}

// copy returns a new copy of m.
pub fn (m Matrix) copy() Matrix {
	mut r := Matrix{
		list: m.list
		dims: m.dims
		offsets: m.offsets
	}
	r.list = []complex.Complex{len: m.list.len}
	arrays.copy[complex.Complex](mut r.list, m.list)
	return r
}

fn (m Matrix) offset(dims []int) int {
	if dims.len != m.dims.len {
		panic('incorrect dimensions')
	}

	mut i := 0
	for n, v in dims {
		if v > m.dims[n] {
			panic('incorrect dimensions')
		}

		i += v * m.offsets[n]
	}

	return i
}

fn (m Matrix) indexes(dims []int) []int {
	mut i := -1
	for n, v in dims {
		if v == -1 {
			if i >= 0 {
				panic('only one dimension index allowed')
			}

			i = n
		} else if v >= m.dims[n] {
			panic('dimension out of bounds')
		}
	}

	if i == -1 {
		panic('must specify one dimension index")')
	}

	mut x := 0
	for n, v in dims {
		if v >= 0 {
			x += m.offsets[n] * v
		}
	}

	mut r := []int{len: m.dims[i]}
	for j, _ in r {
		r[j] = x + m.offsets[i] * j
	}

	return r
}

// dimensions returns the dimension array of the Matrix.
pub fn (m Matrix) dimensions() []int {
	mut r := []int{len: m.dims.len}
	arrays.copy[int](mut r, m.dims)
	return r
}

// dim returns the array of any given index of the Matrix.
// Exactly one value in dims must be -1. This is the array dimension returned.
// For example, using the Matrix documented in new_matrix:
//   m.dim([]int {1, 0, -1}) = []complex.Complex {3, 4, 5, 6}
//   m.dim([]int {0, -1, 2}) = []complex.Complex {3, 7, 1}
//   m.dim([]int {-1, 1, 3}) = []complex.Complex {8, 0}
pub fn (m Matrix) dim(dims []int) []complex.Complex {
	inds := m.indexes(dims)
	mut r := []complex.Complex{len: inds.len}
	for n, v in inds {
		r[n] = m.list[v]
	}

	return r
}

pub fn (mut m Matrix) set_dim(x []complex.Complex, dims []int) {
	inds := m.indexes(dims)
	if x.len != inds.len {
		panic('incorrect array length')
	}

	for n, v in inds {
		m.list[v] = x[n]
	}
}

// value returns the value at the given index.
// m.value([]int {1, 2, 3, 4}) is equivalent to m[1][2][3][4].
pub fn (m Matrix) value(dims []int) complex.Complex {
	return m.list[m.offset(dims)]
}

// set_value sets the value at the given index.
// m.set_value(10, []int {1, 2, 3, 4}) is equivalent to m[1][2][3][4] = 10.
pub fn (mut m Matrix) set_value(x complex.Complex, dims []int) {
	m.list[m.offset(dims)] = x
}

// to_2d returns the 2-D array equivalent of the Matrix.
// Only works on Matrixes of 2 dimensions.
pub fn (m Matrix) to_2d() [][]complex.Complex {
	if m.dims.len != 2 {
		panic('can only convert 2d Matrixes')
	}

	mut r := [][]complex.Complex{len: m.dims[0]}
	for i := 0; i < m.dims[0]; i++ {
		r[i] = []complex.Complex{len: m.dims[1]}
		arrays.copy[complex.Complex](mut r[i], m.list[i * m.dims[1]..(i + 1) * m.dims[1]])
	}

	return r
}

// equal_approx returns true if the Matrixes are very close, else false.
// Comparison done using utils.equal_approx_complex().
pub fn (m Matrix) equal_approx(n Matrix, tol f64) bool {
	for i, v in m.dims {
		if v != n.dims[i] {
			return false
		}
	}

	return equal_approx_complex(m.list, n.list, tol)
}
