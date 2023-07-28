module utils

import math.complex

fn test_new_matrix() {
	tol := 1e-8
	c := [
		// vfmt off
		complex.complex(1, 0),
		complex.complex(2, 0)
		complex.complex(3, 0)
		complex.complex(4, 0)

		complex.complex(5, 0)
		complex.complex(6, 0)
		complex.complex(7, 0)
		complex.complex(8, 0)

		complex.complex(9, 0)
		complex.complex(0, 0)
		complex.complex(1, 0)
		complex.complex(2, 0)


		complex.complex(3, 0)
		complex.complex(4, 0)
		complex.complex(5, 0)
		complex.complex(6, 0)

		complex.complex(7, 0)
		complex.complex(8, 0)
		complex.complex(9, 0)
		complex.complex(0, 0)

		complex.complex(4, 0)
		complex.complex(3, 0)
		complex.complex(2, 0)
		complex.complex(1, 0)
		// vfmt on
	]
	mut m := new_matrix(c, [2, 3, 4])

	check_arr(m.dim([1, 0, -1]), to_complex([f64(3), 4, 5, 6]), tol)
	check_arr(m.dim([0, -1, 2]), to_complex([f64(3), 7, 1]), tol)
	check_arr(m.dim([-1, 1, 3]), to_complex([f64(8), 0]), tol)

	s := to_complex([f64(10), 11, 12])
	i := [1, -1, 3]
	m.set_dim(s, i)
	check_arr(m.dim(i), s, tol)

	v := complex.complex(14, 0)
	m.set_value(v, i)
	check_float(m.value(i), v, tol)
}

fn check_arr(have []complex.Complex, want []complex.Complex, tol f64) {
	assert equal_approx_complex(have, want, tol), 'have: ${have} want: ${want}'
}

fn check_float(have complex.Complex, want complex.Complex, tol f64) {
	assert equal_within_abs_complex(have, want, tol)
}
