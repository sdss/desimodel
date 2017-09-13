# See LICENSE.rst for BSD 3-clause license info
# -*- coding: utf-8 -*-
"""
====================
desimodel.focalplane
====================

Provides utilities to model the DESI focal plane.
"""


import os
import numpy as np
from astropy.io import fits
import astropy.units as u

# Define this here to avoid a problem with Sphinx compilation.
try:
    default_offset = 10.886*u.um
except TypeError:
    default_offset = 10.886


def generate_random_vector_field(rms, exponent, n, seed=None, smoothing=0.02):
    """Generate a pair dx, dy of 2D Gaussian random field.

    The random field is generated with a power spectrum P(k) ~ r ** exponent
    and normalized to the specified RMS value.  Smoothing is applied to minimize
    grid artifacts.

    Parameters
    ----------
    rms : float or astropy quantity
        Desired RMS of the generated field values.
    exponent : float
        Exponent of the power spectrum scaling with radius.
    n : int
        Size of the generated array along each axis.
    seed : int
        Random number seed to use. Generated fields should be portable
        across python versions and platforms.
    smoothing : float
        Length scale for smoothing the generated field expressed
        as a fraction of the full width of the field.  Implemented
        as a Gaussian convolution.  No smoothing is applied when
        smoothing is zero.

    Returns
    -------
    tuple
        Tuple dx, dy of 2D arrays containing the generated Gaussian
        random field values. Arrays will have the same units as the
        rms parameter, if any.
    """
    A = np.zeros((n, n), complex)
    kvec = np.fft.fftfreq(n)
    kx, ky = np.meshgrid(kvec, kvec, sparse=True, copy=False)
    ksq = kx ** 2 + ky ** 2
    m = ksq > 0
    gen = np.random.RandomState(seed=seed)
    phase = 2 * np.pi * gen.uniform(size=(n, n))
    A[m] = (ksq[m] ** (exponent / 2) * gen.normal(size=(n, n))[m] *
            np.exp(1.j * phase[m]))
    if smoothing > 0:
        var = (n * smoothing) ** 2 / 2
        A[m] *= np.exp(-ksq[m] * var) / (2 * np.pi)
    offsets = np.fft.ifft2(A)

    # Rescale to the specified RMS radial offset.
    rescale = rms / np.sqrt(np.var(offsets.real) + np.var(offsets.imag))
    dx = offsets.real * rescale
    dy = offsets.imag * rescale

    return dx, dy


def generate_random_centroid_offsets(rms_offset=default_offset, seed=123):
    """Generate random centroid offsets.

    Calls :func:`generate_random_vector_field` to generate offsets with
    a power spectrum exponent of -1 in the expected format.

    The arrays in the files $DESIMODEL/data/throughput/DESI-0347_random_offset_<n>.fits
    were generated by this method with seeds n=1,2,3.

    The default RMS offset value is the quadrature sum of the achromatic
    contributions in cells B43:B72 of the throughput sheet in DESI-0347-v10.

    Parameters
    ----------
    rms_offset : :class:`astropy.Quantity` instance.
        RMS that the generated offsets should have, including units.
    seed : :class:`int`
        Random number seed to use. Generated offsets should be portable
        across python versions and platforms.

    Returns
    -------
    tuple
        Tuple dx, dy of centroid offset arrays with units.
    """
    return generate_random_vector_field(
        rms_offset, exponent=-1.0, n=256, seed=seed, smoothing=0.02)

_tile_radius_deg = None
_tile_radius_mm = None

def get_tile_radius_mm():
    '''Returns radius in mm to the middle of the outermost positioner'''
    global _tile_radius_mm
    if _tile_radius_mm is None:
        import desimodel.io
        fp = desimodel.io.load_fiberpos()
        p = desimodel.io.load_desiparams()
        _tile_radius_mm = np.sqrt(np.max(fp['X']**2 + fp['Y']**2))  # + p['positioners']['radius_max']

    return _tile_radius_mm

def get_tile_radius_deg():
    '''Returns radius in degrees to the middle of the outermost positioner'''
    global _tile_radius_deg
    if _tile_radius_deg is None:
        import scipy.interpolate
        import desimodel.io
        rmax = get_tile_radius_mm()
        platescale = desimodel.io.load_platescale()
        fn = scipy.interpolate.interp1d(platescale['radius'], platescale['theta'], kind='quadratic')
        _tile_radius_deg = float(fn(rmax))

    return _tile_radius_deg

def get_radius_mm(theta):
    """
    Returns an array of radii in mm given an array of radii in degrees using the platescale data
    relative to the center of the focal plane as (0,0). Supports scalar and vector inputs.
    Parameters
    ----------
    theta: An array that represents the angle from the center of the focal plane
    """
    import scipy.interpolate
    import desimodel.io
    platescale = desimodel.io.load_platescale()
    # Uses a quadratic one-dimensional interpolation to approximate the radius in degrees versus radius in mm
    fn = scipy.interpolate.interp1d(platescale['theta'], platescale['radius'], kind = 'quadratic')
    radius = fn(theta)
    if(np.isscalar(theta)):
        return float(radius)
    else:
        return radius

def get_radius_deg(x, y):
    """
    Returns the radius in degrees given x, y coordinates using the platescale data
    Parameters
    ----------
    x: The x coordinate in mm of a location on the focal plane
    y: The y coordinate in mm of a location on the focal plane
    """
    import scipy.interpolate
    import desimodel.io
    radius = np.sqrt(x**2 + y**2)
    platescale = desimodel.io.load_platescale()
    fn = scipy.interpolate.interp1d(platescale['radius'], platescale['theta'], kind = 'quadratic')
    degree = float(fn(radius))
    return degree

def xy2radec(telra, teldec, x, y):
    """
    Returns the new RA and Dec of an x, y position on the focal plane
    in the sky given an arbitrary telescope pointing in RA and Dec
    Parameters
    ----------
    telra: a float signifying the telescope's RA pointing in degrees
    teldec: a float signifying the telescope's Dec pointing in degrees
    x: The x coordinate in mm of a location on the focal plane
    y: The y coordinate in mm of a location on the focal plane
    """
    from math import atan2, acos, radians, degrees
    
    # radial distance on the focal plane in degrees
    r_deg = get_radius_deg(x, y)
    # q signifies the angle the position makes with the +x-axis of focal plane
    q = np.degrees(np.arctan2(y, x))
    
    coord = np.zeros(shape=(3,1))
    coord[0] = 1
    
    # Clockwise rotation around the z-axis by the radial distance to a point on the focal plane in radians
    zrotate = np.zeros(shape=(3,3))
    r_rad = radians(r_deg)
    zrotate[0] = [np.cos(r_rad), np.sin(r_rad), 0]
    zrotate[1] = [-np.sin(r_rad), np.cos(r_rad), 0]
    zrotate[2] = [0, 0, 1]
    
    # Counter-clockwise rotation around the x-axis
    xrotate = np.zeros(shape=(3,3))
    q_rad = radians(q)
    xrotate[0] = [1, 0, 0]
    xrotate[1] = [0, np.cos(q_rad), -np.sin(q_rad)]
    xrotate[2] = [0, np.sin(q_rad), np.cos(q_rad)]
    
    # Counter-clockwise rotation around y axis by declination of the tile center
    decrotate = np.zeros(shape=(3,3))
    teldec_rad = radians(teldec)
    decrotate[0] = [np.cos(teldec_rad), 0, -np.sin(teldec_rad)]
    decrotate[1] = [0, 1, 0]
    decrotate[2] = [np.sin(teldec_rad), 0, np.cos(teldec_rad)]
    
    # Counter-clockwise rotation around the z-axis by the right ascension of the tile center
    rarotate = np.zeros(shape=(3,3))
    telra_rad = radians(telra)
    rarotate[0] = [np.cos(telra_rad), -np.sin(telra_rad), 0]
    rarotate[1] = [np.sin(telra_rad), np.cos(telra_rad), 0]
    rarotate[2] = [0, 0, 1]
    
    coord1 = np.matmul(zrotate, coord)
    coord2 = np.matmul(xrotate, coord1)
    coord3 = np.matmul(decrotate, coord2)
    coord4 = np.matmul(rarotate, coord3)
    
    ra_rad = atan2(coord4[1], coord4[0])
    dec_rad = (np.pi / 2) - acos(coord4[2] / np.sqrt((coord4[0]**2) + (coord4[1]**2) + (coord4[2]**2)))
    ra_deg = degrees(ra_rad)
    dec_deg = degrees(dec_rad)
    # Value can be 360, which should be 0
    ra = ra_deg % 360
    return ra, dec_deg

def radec2xy(telra, teldec, ra, dec):
    """
    Returns arrays of the x, y positions of given celestial objects
    on the focal plane given an arbitrary telescope pointing in RA and Dec and
    arrays of the RA and Dec of celestial objects in the sky. Implements the Haversine
    formula
    Parameters
    ----------
    telra: a scalar float signifying the telescope's RA pointing in degrees
    teldec: a scalar float signifying the telescope's Dec pointing in degrees
    ra: An array of RA values for locations in the sky
    dec: An array of declination values for locations in the sky
    """
    import numpy as np
    # Inclination is 90 degrees minus the declination in degrees
    dec = np.asarray(dec)
    inc = 90 - dec
    ra = np.asarray(ra)
    #inc = 90 - dec
    x0 = np.sin(np.radians(inc)) * np.cos(np.radians(ra))
    y0 = np.sin(np.radians(inc)) * np.sin(np.radians(ra))
    z0 = np.cos(np.radians(inc))
    coord = [x0, y0, z0]
    
    # Clockwise rotation around y axis by declination of the tile center
    decrotate = np.zeros(shape=(3,3))
    teldec_rad = np.radians(teldec)
    decrotate[0] = [np.cos(teldec_rad), 0, np.sin(teldec_rad)]
    decrotate[1] = [0, 1, 0]
    decrotate[2] = [-np.sin(teldec_rad), 0, np.cos(teldec_rad)]
    
    # Clockwise rotation around the z-axis by the right ascension of the tile center
    rarotate = np.zeros(shape=(3,3))
    telra_rad = np.radians(telra)
    rarotate[0] = [np.cos(telra_rad), np.sin(telra_rad), 0]
    rarotate[1] = [-np.sin(telra_rad), np.cos(telra_rad), 0]
    rarotate[2] = [0, 0, 1]
    
    coord1 = np.matmul(rarotate, coord)
    coord2 = np.matmul(decrotate, coord1)
    x = coord2[0]
    y = coord2[1]
    z = coord2[2]
    
    newteldec = 0
    newtelra = 0
    ra_rad = np.arctan2(y, x)
    dec_rad = (np.pi / 2) - np.arccos(z / np.sqrt((x**2) + (y**2) + (z**2)))
    radius_rad = 2 * np.arcsin(np.sqrt((np.sin((dec_rad - newteldec) / 2)**2) + ((np.cos(newteldec)) * np.cos(dec_rad) * (np.sin((ra_rad - newtelra) / 2)**2))))
    radius_deg = np.degrees(radius_rad)
    
    q_rad = np.arctan2(-z, -y)
    
    radius_mm = get_radius_mm(radius_deg)
    x_focalplane = radius_mm * np.cos(q_rad)
    y_focalplane = radius_mm * np.sin(q_rad)
    
    return x_focalplane, y_focalplane

class FocalPlane(object):
    """A class for modeling the DESI focal plane and converting between
    focal plane coordinates (in mm) and RA, Dec on the sky (in degrees).
    Provides utility functions for mapping which positioners cover
    which (RA, Dec) or (x, y) locations and vice versa.

    Parameters
    ----------
    ra, dec : :class:`float`
        Initialize DESI focal plane model with the telescope pointing
        at (`ra`, `dec`) in degrees.
    """

    def __init__(self, ra=0.0, dec=0.0):
        """
        """
        # Read $DESIMODEL/data/focalplane/fiberpos.fits and platescale.txt
        # to construct focal plane model.  May also need data/desi.yaml .
        self._check_radec(ra, dec)
        self.ra = ra
        self.dec = dec
        self._fiberpos_file = os.path.join(os.environ['DESIMODEL'],
                                           'data', 'focalplane',
                                           'fiberpos.fits')
        with fits.open(self._fiberpos_file) as hdulist:
            self.fiberpos = hdulist[1].data

    def _check_radec(self, ra, dec):
        """Raise ValueError if RA or dec are out of bounds.
        """
        if np.any( (ra < 0) | (ra >= 360) ):
            raise ValueError("RA must be 0 <= RA < 360")
        if np.any( (dec < -90) | (dec > +90) ):
            raise ValueError("Dec must be -90 <= dec <= 90")

    def set_tele_pointing(self, ra, dec):
        """Set telescope pointing to (RA, Dec) in degrees.

        Parameters
        ----------
        ra
        dec : :class:`float`
            Telescope pointing in degrees.
        """
        self._check_radec(ra, dec)
        self.ra = ra
        self.dec = dec

    def plate_dist(self, theta):
        """Returns the radial distance on the plate (mm) given the angle
        (radians). This is a fit to some data, it should be calculated on the
        fly at __init__.

        Parameters
        ----------
        theta : :class:`float`
            Angle in radians.

        Returns
        -------
        :class:`float`
            Radial distance in mm.
        """
        p = np.array([8.297E5, -1750.0, 1.394E4, 0.0])
        radius = 0.0
        for i in range(4):
            radius = theta*radius + p[i]
        return radius

    def plate_angle(self, radius):
        """Returns the angular distance on the plate (radians) given the
        radial distance to the plate (mm).

        It uses a Newton-Raphson method.

        Parameters
        ----------
        radius : :class:`float`
            Plate distance in mm.

        Returns
        -------
        :class:`float`
            Angular distance in radians.
        """
        angle_guess = 0.0 * radius
        dist_guess = self.plate_dist(angle_guess) - radius
        while np.any(np.fabs(dist_guess) > 1E-8):
            derivative = (self.plate_dist(angle_guess+1E-5) -
                          self.plate_dist(angle_guess))/1E-5
            delta_guess = - dist_guess/derivative
            angle_guess = angle_guess + delta_guess
            dist_guess = self.plate_dist(angle_guess) - radius
        return angle_guess

    def radec2xy(self, ra, dec):
        """Convert (RA, Dec) in degrees to (x, y) in mm on the focal plane
        given the current telescope pointing.

        If RA and Dec are floats, returns a tuple (x, y) of floats
        If RA and Dec are numpy arrays, returns a tuple (x, y) of numpy arrays

        Parameters
        ----------
        ra, dec : :class:`float` or :class:`numpy.ndarray`
            Sky position.

        Returns
        -------
        :func:`tuple`
            A tuple containing the (x, y) coordinates in mm.
        """
        self._check_radec(ra, dec)

        object_theta = (90.0 - dec)*np.pi/180.0
        object_phi = ra*np.pi/180.0
        o_hat0 = np.sin(object_theta)*np.cos(object_phi)
        o_hat1 = np.sin(object_theta)*np.sin(object_phi)
        o_hat2 = np.cos(object_theta)

        tile_theta = (90.0 - self.dec)* np.pi/180.0
        tile_phi = self.ra * np.pi/180.0
        t_hat0 = np.sin(tile_theta)*np.cos(tile_phi)
        t_hat1 = np.sin(tile_theta)*np.sin(tile_phi)
        t_hat2 = np.cos(tile_theta)

        # We make a rotation on o_hat, so that t_hat ends up aligned with
        # the unit vector along z. This is composed by a first rotation around
        # z of an angle pi/2 - phi and a second rotation around x by an angle
        # theta, where theta and phi are the angles describin t_hat.

        costheta = t_hat2
        sintheta = np.sqrt(1.0-costheta*costheta) + 1E-12
        cosphi = t_hat0/sintheta
        sinphi = t_hat1/sintheta
        # First rotation, taking into account that cos(pi/2 -phi) =
        # sin(phi) and sin(pi/2-phi)=cos(phi)
        n_hat0 = sinphi*o_hat0 - cosphi*o_hat1
        n_hat1 = cosphi*o_hat0 + sinphi*o_hat1
        n_hat2 = o_hat2
        # Second rotation
        nn_hat0 = n_hat0
        nn_hat1 = costheta*n_hat1 - sintheta*n_hat2
        nn_hat2 = sintheta*n_hat1 + costheta*n_hat2
        # Now find the radius on the plate
        theta = np.sqrt(nn_hat0*nn_hat0 + nn_hat1*nn_hat1)
        radius = self.plate_dist(theta)
        x = radius * nn_hat0/theta
        y = radius * nn_hat1/theta
        return (x, y)

    def xy2radec(self, x, y):
        """Convert (x, y) in mm on the focal plane to (ra_object, dec_object)
        in degrees on the sky given the current telescope pointing towards
        (RA, Dec).

        `x`, `y` must be floats. This function is vectorized in xy2radec(),
        which doesn't appear to exist.

        Parameters
        ----------
        x, y : :class:`float`
            Position on the focal plane in mm.

        Returns
        -------
        :func:`tuple`
            Coordinates of object.
        """
        self._check_radec(self.ra, self.dec)
        # The implementation of this function takes the same conventions
        # as radec2xy to make rotations, only that it happens in the reverse
        # order.
        #
        # This is the final position of the tile vector, which starts
        # parallel to  z_hat.
        tile_theta = (90.0 - self.dec)*np.pi/180.0
        tile_phi = self.ra*np.pi/180.0
        t_hat0 = np.sin(tile_theta)*np.cos(tile_phi)
        t_hat1 = np.sin(tile_theta)*np.sin(tile_phi)
        t_hat2 = np.cos(tile_theta)
        # Define sin and cos of the angles for the final tile vector.
        costheta = t_hat2
        sintheta = np.sqrt(1.0-costheta*costheta) + 1E-10
        cosphi = t_hat0/sintheta
        sinphi = t_hat1/sintheta
        # Find the initial position of the object vector when the tile
        # starts parallel to z_hat.
        radius = np.sqrt(x*x + y*y)
        object_theta = self.plate_angle(radius)
        object_phi = np.arctan2(y, x)
        o_hat0 = np.sin(object_theta)*np.cos(object_phi)
        o_hat1 = np.sin(object_theta)*np.sin(object_phi)
        o_hat2 = np.cos(object_theta)
        # First rotation, around x by an angle -theta.
        n_hat0 = o_hat0
        n_hat1 = costheta*o_hat1 + sintheta*o_hat2
        n_hat2 = -sintheta*o_hat1 + costheta*o_hat2
        # Second rotation around z_hat by -(pi/2-phi), taking into
        # account that cos(pi/2 -phi) = sin(phi) and
        # sin(pi/2-phi)=cos(phi).
        nn_hat0 = sinphi*n_hat0 + cosphi*n_hat1
        nn_hat1 = -cosphi*n_hat0 + sinphi*n_hat1
        nn_hat2 = n_hat2
        # Convert from unit vectors to ra, dec.
        object_theta = np.arccos(nn_hat2)
        object_phi = np.arctan2(nn_hat1, nn_hat0)
        object_dec = 90.0 - (180.0/np.pi)*object_theta
        # Due to rounding imprecisions the remainder has to be taken
        # two times (!).
        object_ra = np.remainder((object_phi*180.0/np.pi), 360.0)
        object_ra = np.remainder(object_ra, 360.0)
        self._check_radec(object_ra, object_dec)
        return (object_ra, object_dec)

    def radec2pos(self, ra, dec):
        """Identify which positioners cover (`ra`, `dec`).

        If `ra`, `dec` are floats, return an array of positioner IDs that
        cover it. The array could be empty if no positioner covers that
        location.

        If `ra`, `dec` are numpy arrays, return a list of arrays.
        The ith element is an array of positioner IDs that cover
        (ra[i], dec[i]).

        .. warning:: This method is not implemented!
        """
        raise NotImplementedError


# Not sure what to call it, but it would also be useful to have a
# function that takes (ra, dec) arrays and returns a list with one
# element per positioner, giving the indices of the (ra,dec) arrays
# that that positioner covers.
#
# ditto for xy2pos(), xy2fiber(), etc.
# positioner = thing on the focal plane
# fiber = numerically increasing in order on the spectrograph CCDs

def fiber_area_arcsec2(x, y):
    '''
    Returns area of fibers at (x,y) in arcsec^2
    '''
    from desimodel.io import load_desiparams, load_platescale
    params = load_desiparams()
    fiber_dia = params['fibers']['diameter_um']
    x = np.asarray(x)
    y = np.asarray(y)
    r = np.sqrt(x**2 + y**2)

    #- Platescales in um/arcsec
    ps = load_platescale()
    radial_scale = np.interp(r, ps['radius'], ps['radial_platescale'])
    az_scale = np.interp(r, ps['radius'], ps['az_platescale'])

    #- radial and azimuthal fiber radii in arcsec
    rr  = 0.5 * fiber_dia / radial_scale
    raz = 0.5 * fiber_dia / az_scale
    fiber_area = (np.pi * rr * raz)
    return fiber_area

def on_gfa(telra, teldec, ra, dec, buffer_arcsec = 100):
    """
    Checks if a target is on any of the 10 GFAs given telra, teldec and an array of RA and Dec pointings,
    as well as a parameter for degrees of tolerance one would like to allow. When using
    desimodel.footprint.find_points_in_tiles(tiles, ra, dec, radius) with this function to
    check what points are on the GFAs, the default radius parameter should be set to 1.65 (degrees),
    so that boundary GFA area actually encompasses points normally outside of the tile.
    Parameters:
    telra: The telescope's arbitrary RA pointing
    teldec: The telescope's arbitrary Dec pointing
    ra: An array of RA values for locations in the sky
    dec: An array of declination values for locations in the sky
    buffer_arcsec: A value in arcseconds on the sky of how much tolerance
    one would allow for seeing if a target is on the gfa.
    Returns:
    targetindices: a list  of targets with their respective indices in the
    RA and Dec list passed in that fall on certain GFAs denoted by the index
    in the gfaindices list.
    gfaindices: a list equal in length with the targetindices list with the gfa location 0-9 as each element
    """
    import desimodel.footprint
    # If any calculated area is under the threshold area, it is mathematically impossible
    THRESHOLD_AREA = 469.7
    MIN_TOLERANCE = 0.001
    
    inrangeindices = desimodel.footprint.find_points_radec(telra, teldec, ra, dec, 1.65)
    if not inrangeindices:
        return np.array([]), np.array([])
    inrangeindices = np.asarray(inrangeindices)
    
    targetx, targety = desimodel.focalplane.radec2xy(telra, teldec, ra[inrangeindices], dec[inrangeindices])
    
    x_tolerance, y_tolerance = _degrees2xytolerance(buffer_arcsec)
    
    targetindices = []
    gfaindices = []
    
    # x and y hold the 40 new GFA coordinates
    x, y = _shift_gfa_points(x_tolerance, y_tolerance)
    # The area boundary's value is the area of the gfa plus some tolerance.
    AREA_BOUNDARY = _retrieve_minimum_boundary(x_tolerance, y_tolerance) + MIN_TOLERANCE

    targetx = np.asarray(targetx)
    targety = np.asarray(targety)
    # Method to check if point is inside the rectangle
    for gfaid in range(0, 40, 4):
        # a1 through a4 are edge lengths of the rectangle formed by corners of the GFAs
        a1 = np.sqrt((x[gfaid] - x[gfaid + 1])**2 + (y[gfaid] - y[gfaid + 1])**2)
        a2 = np.sqrt((x[gfaid + 1] - x[gfaid + 2])**2 + (y[gfaid + 1] - y[gfaid + 2])**2)
        a3 = np.sqrt((x[gfaid + 2] - x[gfaid + 3])**2 + (y[gfaid + 2] - y[gfaid + 3])**2)
        a4 = np.sqrt((x[gfaid + 3] - x[gfaid])**2 + (y[gfaid + 3] - y[gfaid])**2)
        # b1 through b4 are the line segments from each corner to the target location
        b1 = np.sqrt((x[gfaid] - targetx)**2 + (y[gfaid] - targety)**2)
        b2 = np.sqrt((x[gfaid + 1] - targetx)**2 + (y[gfaid + 1] - targety)**2)
        b3 = np.sqrt((x[gfaid + 2] - targetx)**2 + (y[gfaid + 2] - targety)**2)
        b4 = np.sqrt((x[gfaid + 3] - targetx)**2 + (y[gfaid + 3] - targety)**2)
        # Calculating areas of triangles using Heron's Formula
        u1 = (a1 + b1 + b2) / 2.0
        u2 = (a2 + b2 + b3) / 2.0
        u3 = (a3 + b3 + b4) / 2.0
        u4 = (a4 + b4 + b1) / 2.0
        area1 = np.sqrt((u1 * (u1 - a1) * (u1 - b1) * (u1 - b2)).clip(0))
        area2 = np.sqrt((u2 * (u2 - a2) * (u2 - b2) * (u2 - b3)).clip(0))
        area3 = np.sqrt((u3 * (u3 - a3) * (u3 - b3) * (u3 - b4)).clip(0))
        area4 = np.sqrt((u4 * (u4 - a4) * (u4 - b4) * (u4 - b1)).clip(0))
        targetarea = area1 + area2 + area3 + area4
         
                     
        assert np.all(targetarea > THRESHOLD_AREA)
                     
        if(any(targetarea < AREA_BOUNDARY) and all(targetarea > THRESHOLD_AREA)):
            newtargetindices = np.where(targetarea < AREA_BOUNDARY)
            targetindices.extend(newtargetindices[0])
            gfaindices.extend([int(gfaid / 4)] * len(newtargetindices[0]))
    return inrangeindices[targetindices], gfaindices

def _retrieve_minimum_boundary(x_tolerance, y_tolerance):
    """
    Used as a helper function to the on_gfa function to find the minimum boundary
    area for a point to lie inside a certain GFA given an tolerance in x and y in mm
    Parameters:
    x_tolerance: tolerance in x in mm
    y_tolerance: tolerance in y in mm
    Returns:
    targetarea: the minimum boundary area for the procedure to check if a point is inside the GFA
    """
    import desimodel.footprint
    import desimodel.focalplane
    
    targetx = 116.279135121
    targety = -372.885546514
    #6.644525362152656, -9.055425745149217 GUARANTEED TO BE IN GFA (RA, DEC)
    #x, y = desimodel.focalplane.radec2xy(7.11, -10.53, targetx, targety)
    # If any calculated area is under the threshold area, it is mathematically impossible
    THRESHOLD_AREA = 469.7
    MIN_TOLERANCE = 0.001
    # The area boundary's value is the area of the gfa plus some tolerance.
    
    # x and y hold the 40 new GFA coordinates
    x, y = _shift_gfa_points(x_tolerance, y_tolerance)
    
    targetx = np.asarray(targetx)
    targety = np.asarray(targety)
    # Method to check if point is inside the rectangle
    for gfaid in range(0, 4, 4):
        # a1 through a4 are edge lengths of the rectangle formed by corners of the GFAs
        a1 = np.sqrt((x[gfaid] - x[gfaid + 1])**2 + (y[gfaid] - y[gfaid + 1])**2)
        a2 = np.sqrt((x[gfaid + 1] - x[gfaid + 2])**2 + (y[gfaid + 1] - y[gfaid + 2])**2)
        a3 = np.sqrt((x[gfaid + 2] - x[gfaid + 3])**2 + (y[gfaid + 2] - y[gfaid + 3])**2)
        a4 = np.sqrt((x[gfaid + 3] - x[gfaid])**2 + (y[gfaid + 3] - y[gfaid])**2)
        # b1 through b4 are the line segments from each corner to the target location
        b1 = np.sqrt((x[gfaid] - targetx)**2 + (y[gfaid] - targety)**2)
        b2 = np.sqrt((x[gfaid + 1] - targetx)**2 + (y[gfaid + 1] - targety)**2)
        b3 = np.sqrt((x[gfaid + 2] - targetx)**2 + (y[gfaid + 2] - targety)**2)
        b4 = np.sqrt((x[gfaid + 3] - targetx)**2 + (y[gfaid + 3] - targety)**2)
        # Calculating areas of triangles using Heron's Formula
        u1 = (a1 + b1 + b2) / 2.0
        u2 = (a2 + b2 + b3) / 2.0
        u3 = (a3 + b3 + b4) / 2.0
        u4 = (a4 + b4 + b1) / 2.0
        area1 = np.sqrt(u1 * (u1 - a1) * (u1 - b1) * (u1 - b2))
        area2 = np.sqrt(u2 * (u2 - a2) * (u2 - b2) * (u2 - b3))
        area3 = np.sqrt(u3 * (u3 - a3) * (u3 - b3) * (u3 - b4))
        area4 = np.sqrt(u4 * (u4 - a4) * (u4 - b4) * (u4 - b1))
        targetarea = area1 + area2 + area3 + area4
         
        assert np.all(targetarea > THRESHOLD_AREA)
        return targetarea

def _degrees2xytolerance(buffer_arcsec):
    """
    Used as a helper function to the on_gfa function to find the tolerance in x and y
    given a tolerance in arcseconds
    Parameters:
    buffer_arcsec: a tolerance in arcseconds for checking if a point is on the GFA
    Returns:
    x_tolerance: tolerance in x in mm
    y_tolerance: tolerance in y in mm
    """
    # Uses the center of a given GFA from DESI-0530-v13 Excel Spreadsheet to find the tolerance
    import desimodel.io
    import scipy.interpolate
    platescale = desimodel.io.load_platescale()
    fn = scipy.interpolate.interp1d(platescale['radius'], platescale['radial_platescale'], kind = 'quadratic')
    fn1 = scipy.interpolate.interp1d(platescale['radius'], platescale['az_platescale'], kind = 'quadratic')
    # Center of a given GFA from DESI-0530-v13 Excel Spreadsheet
    x = 333.738
    y = 217.766
    radius = np.sqrt(x**2 + y**2)
    # Platescales are in units of microns per arcsecond
    r_ps = fn(radius)
    az_ps = fn(radius)
    x_tolerance = buffer_arcsec / (10**3) * r_ps
    y_tolerance = buffer_arcsec / (10**3) * az_ps
    return x_tolerance, y_tolerance

def _shift_gfa_points(deltax, deltay):
    """
    Used as a helper function to the on_gfa function to find the new
    GFA locations after incorporating a tolerance in x and y
    Parameters:
    deltax: tolerance in x in mm
    deltay: tolerance in y in mm
    Returns:
    Returns the 40 new GFA locations in x and y
    """
    import numpy as np
    x = [-125.10482863, -129.83038525, -159.04283509, -154.31646944]
    y = [-370.01790486, -384.56223777, -375.05643893, -360.51151824]
    point1 = [x[2], y[2]]
    point2 = [x[1], y[1]]
    vector1 = [(point2[0] - point1[0]), (point2[1] - point1[1])]
    vector2 = [1, 0]
    # Angle between vector1 and vector 2 using dot product
    angle = np.arccos((np.dot(vector1, vector2))/(np.sqrt((vector1[0]**2) + (vector1[1]**2))))
    
    shiftmat = np.zeros(shape=(2,2))
    shiftmat[0] = [np.cos(angle), -np.sin(angle)]
    shiftmat[1] = [np.sin(angle), np.cos(angle)]
    reverseshift= np.zeros(shape=(2,2))
    reverseshift[0] = [np.cos(angle), np.sin(angle)]
    reverseshift[1] = [-np.sin(angle), np.cos(angle)]
    
    # Shifts the initial coordinates to be parallel to the vector [1, 0]
    coord = np.zeros(shape=(2,1))
    oldxcoord = x
    oldycoord = y
    for i in range(4):
        coord[0] = oldxcoord[i]
        coord[1] = oldycoord[i]
        newcoord = np.matmul(shiftmat, coord)
        oldxcoord[i] = newcoord[0]
        oldycoord[i] = newcoord[1]
        if(i == 0 or i == 1):
            x[i] = newcoord[0] + deltax
        else:
            x[i] = newcoord[0] - deltax
        if(i == 1 or i == 2):
            y[i] = newcoord[1] - deltay
        else:
            y[i] = newcoord[1] + deltay
    oldxcoord = x
    oldycoord = y
    for i in range(4):
        coord[0] = oldxcoord[i]
        coord[1] = oldycoord[i]
        newcoord = np.matmul(reverseshift, coord)
        oldxcoord[i] = newcoord[0]
        oldycoord[i] = newcoord[1]
        x[i] = newcoord[0]
        y[i] = newcoord[1]
    
    rotatemat = np.zeros(shape=(2,2))
    rotatemat[0] = [np.cos(np.radians(36)), -np.sin(np.radians(36))]
    rotatemat[1] = [np.sin(np.radians(36)), np.cos(np.radians(36))]
    return _find_new_gfa_coordinates(x, y, rotatemat)

def _find_new_gfa_coordinates(x, y, rotatemat):
    """
    Used as a helper function to the on_gfa function to find the new
    GFA coordinates given a list of x coordinates, y coordinates, and a rotation matrix
    Parameters:
    x: a list of x coordinates for the GFAs
    y: a list of y coordinates for the GFAs
    rotatemat: a matrix for rotating the respective coordinates
    Returns:
    x_all: a complete list of the 40 GFA x coordinates
    y_all: a complete list of the 40 GFA y coordinates
    """
    import numpy as np
    x_all = np.zeros(shape=(40,1))
    y_all = np.zeros(shape=(40,1))
    coord = np.zeros(shape=(2,1))
    gfacoord = np.zeros(shape=(4, 2))
    oldxcoord = x
    oldycoord = y
    counter = 0
    for j in range(10):
        for i in range(4):
            coord[0] = oldxcoord[i]
            coord[1] = oldycoord[i]
            newcoord = np.matmul(rotatemat, coord)
            oldxcoord[i] = newcoord[0]
            oldycoord[i] = newcoord[1]
            gfacoord[i] = [newcoord[0], newcoord[1]]
            x_all[counter] = newcoord[0]
            y_all[counter] = newcoord[1]
            counter += 1
    return x_all, y_all

def on_tile_gfa(tileid, targets, buffer_arcsec = 100):
    """
    This function takes a tileid, a table of targets, and an optional
    buffer_arcsec parameter to return the indices of targets lying on the GFA
    as well as the GFA locations from 0-9
    Parameters:
    tileid: (int) DESI tile ID, used to lookup telescope (RA, dec)
    targets: table with columns RA, DEC
    Options:
    buffer_arcsec: (float) additional buffer region around GFA to include
    Returns:
    targetindices: list of indices for targets that are covered by GFA number
    in corresponding gfaindices
    gfaindices: list of indices corresponding to 0-9 GFA location
    """
    import desimodel.footprint
    telra, teldec = desimodel.footprint.get_tile_radec(tileid)
    return on_gfa(telra, teldec, targets['RA'], targets['DEC'], buffer_arcsec)

def get_gfa_targets(targets, rfluxlim = 1000, tiles = None, buffer_arcsec = 100):
    """
    This function takes a table of targets, as well as optional parameters
    including a minimum flux in the r-band, a list of tiles, and a buffer in arcseconds
    and returns a table of targets on the GFA satisfying a minimum flux_r
    Parameters:
    targets: table with columns RA, DEC, FLUX_R
    Options:
    rfluxlim: (float) r-band flux limit; default 1000 = rmag 15
    tiles: table of tiles, default to desimodel.io.load_tiles()
    buffer_arcsec: (float) additional buffer region around GFA to include
    Returns subset of input `targets` with additional columns:
    TILEID: (integer) DESI tile ID
    GFA_LOC: (integer) GFA location [0-9]
    Note that the same target could be repeated with different TILEID, GFA_LOC
    Note also that the function returns an empty list if no targets are on any GFAs or of sufficient brightness
    """
    if(tiles is None):
        import desimodel.io
        tiles = desimodel.io.load_tiles()
    import desimodel.footprint
    points = desimodel.footprint.find_points_in_tiles(tiles, targets['RA'], targets['DEC'])
    alltargetindices = []
    tileidlist = []
    gfaidlist = []
    # Checks if the flux_r meets a minimum threshold
    brightindices = np.where(targets['FLUX_R'] > rfluxlim)
    if(brightindices[0].size == 0):
        return []
    counter = 0
    for lists in points:
        if lists:
            tileid = tiles[counter]['TILEID']
            targetindices, gfaindices = on_tile_gfa(tileid, targets[brightindices[0]], buffer_arcsec)
            tileidlist.extend([tileid] * len(targetindices))
            alltargetindices.extend(targetindices)
            gfaidlist.extend(gfaindices)
        counter += 1
    validtargets = targets[brightindices[0]][alltargetindices]
    tileidlist = np.asarray(tileidlist)
    gfaidlist = np.asarray(gfaidlist)
    validtargets['TILEID'] = tileidlist
    validtargets['GFA_LOC'] = gfaidlist
    return validtargets

