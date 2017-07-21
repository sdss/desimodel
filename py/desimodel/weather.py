# See LICENSE.rst for BSD 3-clause license info
# -*- coding: utf-8 -*-
"""
=================
desimodel.weather
=================

Model of the expected weather conditions at KPNO during the DESI survey.

See DESI-doc-3087 for details.
"""
from __future__ import print_function, division

import numpy as np

import scipy.interpolate
import scipy.special


def whiten_transforms_from_cdf(x, cdf):
    """
    Calculate a pair of transforms to whiten and unwhiten a distribution.

    The whitening transform is monotonic and invertible.

    Parameters
    ----------
    x : array
        1D array of non-decreasing values giving bin edges for the distribution
        to whiten and unwhiten.
    cdf : array
        1D array of non-decreasing values giving the cummulative probability
        density associated with each bin edge.  Does not need to be normalized.
        Must have the same length as x.

    Returns
    -------
    tuple
        Tuple (F,G) of callable objects that whiten y=F(x) and unwhiten x=G(y)
        samples x of the input distribution, so that y has a Gaussian
        distribution with zero mean and unit variance.
    """
    x = np.asarray(x)
    cdf = np.asarray(cdf)
    if x.shape != cdf.shape:
        raise ValueError('Input arrays must have same shape.')
    if len(x.shape) != 1:
        raise ValueError('Input arrays must be 1D.')
    if not np.all(np.diff(x) >= 0):
        raise ValueError('Values of x must be increasing.')
    if not np.all(np.diff(cdf) >= 0):
        raise ValueError('Values of cdf must be increasing.')
    # Normalize.
    cdf /= cdf[-1]
    # Use linear interpolation for the forward and inverse transforms between
    # the input range and Gaussian CDF values.
    args = dict(
        kind='linear', assume_sorted=True, copy=False, bounds_error=True)
    forward = scipy.interpolate.interp1d(x, cdf, **args)
    backward = scipy.interpolate.interp1d(cdf, x, **args)
    # Add wrappers to convert between CDF and PDF samples.
    root2 = np.sqrt(2)
    forward_transform = (
        lambda x: root2 * scipy.special.erfinv(2 * forward(x) - 1))
    inverse_transform = (
        lambda y: backward(0.5 * (1 + scipy.special.erf(y / root2))))
    return forward_transform, inverse_transform


def whiten_transforms(data, data_min=None, data_max=None):
    """Calculate a pair of transforms to whiten and unwhiten a distribution.

    Uses :func:`desimodel.weather.whiten_transforms_from_cdf`.

    Parameters
    ----------
    data : array
        1D array of samples from the distribution to whiten.
    data_min : float or None
        Clip the distribution to this minimum value, or at min(data) if None.
        Must be <= min(data).
    data_max : float or None
        Clip the distribution to this maximum value, or at max(data) if None.
        Must be >= max(data).

    Returns
    -------
    tuple
        See :func:`desimodel.weather.whiten_transforms_from_cdf`.
    """
    n_data = len(data)
    # Sort the input data with padding at each end for the min/max values.
    sorted_data = np.empty(shape=n_data + 2, dtype=data.dtype)
    sorted_data[1:-1] = np.sort(data)
    if data_min is None:
        sorted_data[0] = sorted_data[1]
    else:
        if data_min > sorted_data[1]:
            raise ValueError('data_min > min(data)')
        sorted_data[0] = data_min
    if data_max is None:
        sorted_data[-1] = sorted_data[-2]
    else:
        if data_max < sorted_data[-2]:
            raise ValueError('data_max < max(data)')
        sorted_data[-1] = data_max
    # Calculate the Gaussian CDF value associated with each input value in
    # sorted order. The pad values are associated with CDF = 0, 1 respectively.
    cdf = np.arange(n_data + 2) / (n_data + 1.)

    return whiten_transforms_from_cdf(sorted_data, cdf)


def _seeing_fit_model(x):
    """Evalute the fit to MzLS seeing described in DESI-doc-3087.
    """
    p = np.array([  0.07511146,   0.44276671,  23.02442192,  38.07691498])
    y = (1 + ((x - p[0]) / p[1]) ** 2) ** (-p[2]) * x ** p[3]
    return y / (y.sum() * np.gradient(x))


def get_seeing_pdf(median_seeing=1.1, max_seeing=2.5, n=250):
    """Return PDF of FWHM seeing for specified clipped median value.

    Note that this is atmospheric seeing, not delivered image quality.
    See DESI-doc-3087 for details.

    Scales the clipped MzLS seeing PDF in order to achieve the requested
    median value.  Note that clipping is applied before scaling, so
    the output PDF is clipped at scale * max_seeing.

    Parameters
    ----------
    median_seeing : float
        Target FWHM seeing value in arcsec. Must be in the range [0.95, 1.30].
    max_seeing : float
        Calculate scaled median using unscaled values below this value.

    Returns
    -------
    tuple
        Tuple (fwhm, pdf) that tabulates pdf[fwhm]. Normalized so that
        np.sum(pdf * np.gradient(fwhm)) = 1.
    """
    # Tabulate the nominal (scale=1) seeing PDF.
    fwhm = np.linspace(0., max_seeing, n)
    pdf = _seeing_fit_model(fwhm)
    pdf /= (pdf.sum() * np.gradient(fwhm))
    cdf = np.cumsum(pdf)
    cdf /= cdf[-1]
    # Tabulate the median as a function of FWHM scale.
    scale = np.linspace(0.9, 1.4, 11)
    median = np.empty_like(scale)
    for i, s in enumerate(scale):
        median[i] = np.interp(0.5, cdf, s * fwhm)
    if median_seeing < median[0] or median_seeing > median[-1]:
        raise ValueError('Requested median is outside allowed range.')
    # Interpolate to find the scale factor that gives the requested median.
    s = np.interp(median_seeing, median, scale)
    return fwhm * s, pdf / s
