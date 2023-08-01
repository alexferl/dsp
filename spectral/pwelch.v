module spectral

import math
import math.complex
import fft
import utils
import window

pub struct PwelchOptions {
	// nfft is the number of data points used in each block for the FFT. Must be
	// even; a power of 2 is most efficient. This should *NOT* be used to get zero
	// padding, or the scaling of the result will be incorrect. Use pad for
	// this instead.
	//
	// The default value is 256.
	nfft int = 256
	// window is a function that returns an array of window values the length
	// of its input parameter. Each segment is scaled by these values.
	//
	// The default is window.Hann, from the dsp/window package.
	window fn (int) []f64 = window.hann
	// pad is the number of points to which the data segment is padded when
	// performing the FFT. This can be different from nfft, which specifies the
	// number of data points used. While not increasing the actual resolution of
	// the psd (the minimum distance between resolvable peaks), this can give
	// more points in the plot, allowing for more detail.
	//
	// The value default is 0, which sets pad equal to nfft.
	pad int
	// noverlap is the number of points of overlap between blocks.
	//
	// The default value is 0 (no overlap).
	noverlap int
	// scale_off specifies whether the resulting density values should be scaled by the
	// scaling frequency, which gives density in units of Hz^-1. This allows for
	// integration over the returned frequency values. The default is set for
	// MATLAB compatibility. Note that this is the opposite of matplotlib style,
	// but with equivalent defaults.
	//
	// The default value is false (enable scaling).
	scale_off bool
}

// pwelch estimates the power spectral density of x using Welch's method.
// fs is the sampling frequency (samples per time unit) of x. fs is used
// to calculate frequencies.
// Returns the power spectral density pxx and corresponding frequencies.
// Designed to be similar to the matplotlib implementation below.
// Reference: http://matplotlib.org/api/mlab_api.html#matplotlib.mlab.psd
// See also: http://www.mathworks.com/help/signal/ref/pwelch.html
pub fn pwelch(mut x []f64, fs f64, o PwelchOptions) ([]f64, []f64) {
	if x.len == 0 {
		return []f64{}, []f64{}
	}

	nfft := o.nfft
	mut pad := o.pad
	noverlap := o.noverlap
	wf := o.window
	enable_scaling := !o.scale_off

	if pad == 0 {
		pad = nfft
	}

	if x.len < nfft {
		x = utils.zero_pad[f64](x, nfft)
	}

	lp := pad / 2 + 1
	scale := f64(2)

	mut segs := segment(x, nfft, noverlap)

	mut pxx := []f64{len: lp}
	for mut seg in segs {
		seg = utils.zero_pad[f64](seg, pad)
		window.apply(mut seg, wf)

		pgram := fft.fft_real(seg)

		for j, _ in pxx {
			mut d := (conj((pgram[j])) * pgram[j]).re / f64(segs.len)

			if j > 0 && j < lp - 1 {
				d *= scale
			}

			pxx[j] += d
		}
	}

	w := wf(nfft)
	mut norm := f64(0)
	for i in w {
		norm += math.pow(i, 2)
	}

	if enable_scaling {
		norm *= fs
	}

	for i, _ in pxx {
		pxx[i] /= norm
	}

	mut freqs := []f64{len: lp}
	coef := fs / f64(pad)
	for i, _ in freqs {
		freqs[i] = f64(i) * coef
	}

	return pxx, freqs
}

fn conj(x complex.Complex) complex.Complex {
	return complex.complex(x.re, -x.im)
}
