# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
lvmmodel.io
============

I/O utility functions for files in lvmmodel.
"""
import os
from astropy.io import fits
import yaml
import numpy as np
import warnings

from lvmutil.log import get_logger
log = get_logger()


_thru = dict()


def load_throughput(channel):
    """Returns specter Throughput object for the given channel 'b', 'r', or 'z'.

    Parameters
    ----------
    channel : {'b', 'r', 'z'}
        Spectrograph channel.
    """
    import specter.throughput
    channel = channel.lower()
    global _thru
    if channel not in _thru:
        thrufile = os.path.join(os.environ['LVMMODEL'], 'data', 'throughput', 'thru-{0}.fits'.format(channel))
        _thru[channel] = specter.throughput.load_throughput(thrufile)
    return _thru[channel]


_psf = dict()


def load_psf(channel):
    """Returns specter PSF object for the given channel 'b', 'r', or 'z'.

    Parameters
    ----------
    channel : {'b', 'r', 'z'}
        Spectrograph channel.
    """
    import specter.psf
    channel = channel.lower()
    global _psf
    if channel not in _psf:
        psffile = os.path.join(os.environ['LVMMODEL'], 'data', 'specpsf', 'psf-{0}.fits'.format(channel))
        _psf[channel] = specter.psf.load_psf(psffile)
    return _psf[channel]


_params = None


def load_lvmparams(config='lvm', telescope='1m'):
    """Returns LVM parameter dictionary loaded from lvmmodel/data/lvm.yaml.

    Parameters:
        config (str):
            Which config yaml to load

        telescope (str):
            Which telescope config to load.
    """

    # build param name
    if config == 'lvm':
        config_name = '{0}_{1}.yaml'.format(config, telescope)
    else:
        config_name = '{0}.yaml'.format(config)

    global _params
    sametele = _params is not None and 'telescope' in _params and telescope == _params['telescope']

    if _params is None or not sametele:
        lvmparamsfile = os.path.join(os.environ['LVMMODEL'], 'data', config_name)
        with open(lvmparamsfile) as par:
            _params = yaml.load(par)

        # - add config and telescope name
        _params['config_name'] = config_name
        _params['telescope'] = telescope

        # - for temporary backwards compability after 'exptime' -> 'exptime_dark'
        if ('exptime' not in _params) and ('exptime_dark' in _params):
            _params['exptime'] = _params['exptime_dark']

        # - Augment params with wavelength coverage from specpsf files
        # - wavemin/max = min/max wavelength covered by *any* fiber on the CCD
        # - wavemin/max_all = min/max wavelength covered by *all* fibers
        for channel in ['b', 'r', 'z']:
            hdr = fits.getheader(findfile('specpsf/psf-{}.fits'.format(channel)), 0)
            _params['ccd'][channel]['wavemin'] = hdr['WAVEMIN']
            _params['ccd'][channel]['wavemax'] = hdr['WAVEMAX']
            _params['ccd'][channel]['wavemin_all'] = hdr['WMIN_ALL']
            _params['ccd'][channel]['wavemax_all'] = hdr['WMAX_ALL']

    return _params


# Added and still needs to be committed and pushed to desihub
_gfa = None


def load_gfa():
    """Returns GFA table from lvmmodel/data/focalplane/gfa.ecsv"""
    global _gfa
    from astropy.table import Table
    # os is imported already in the lvmmodel io.py
    import os
    if _gfa is None:
        gfaFile = os.path.join(os.environ['LVMMODEL'], 'data', 'focalplane', 'gfa.ecsv')
        _gfa = Table.read(gfaFile, format='ascii.ecsv')
    return _gfa


_fiberpos = None


def load_fiberpos():
    """Returns fiberpos table from lvmmodel/data/focalplane/fiberpos.fits.
    """
    global _fiberpos
    from astropy.table import Table
    if _fiberpos is None:
        fiberposfile = os.path.join(os.environ['LVMMODEL'], 'data', 'focalplane', 'fiberpos.fits')
        _fiberpos = Table.read(fiberposfile)
        # - Convert to upper case if needed
        # - Make copy of colnames b/c they are updated during iteration
        for col in list(_fiberpos.colnames):
            if col.islower():
                _fiberpos.rename_column(col, col.upper())

        # - Temporary backwards compatibility for renamed columns
        if 'POSITIONER' in _fiberpos.colnames:
            import warnings
            warnings.warn('old fiberpos.fits with POSITIONER column instead of LOCATION; please update your $LVMMODEL checkout', DeprecationWarning)
            _fiberpos['LOCATION'] = _fiberpos['POSITIONER']
        else:
            _fiberpos['POSITIONER'] = _fiberpos['LOCATION']

        if 'SPECTROGRAPH' in _fiberpos.colnames:
            import warnings
            warnings.warn('old fiberpos.fits with SPECTROGRAPH column instead of SPECTRO; please update your $LVMMODEL checkout', DeprecationWarning)
            _fiberpos['SPECTRO'] = _fiberpos['SPECTROGRAPH']
        else:
            _fiberpos['SPECTROGRAPH'] = _fiberpos['SPECTRO']

    return _fiberpos


_tiles = dict()


def load_tiles(onlydesi=True, extra=False, tilesfile=None, cache=True):
    """Return DESI tiles structure from lvmmodel/data/footprint/desi-tiles.fits.

    Parameters
    ----------
    onlydesi : :class:`bool` (default True)
        If ``True``, trim to just the tiles in the DESI footprint.
    extra : :class:`bool`, (default False)
        If ``True``, include extra layers with PROGRAM='EXTRA'.
    tilesfile : (str)
        Name of tiles file to load; or None for default.
        Without path, look in $LVMMODEL/data/footprint, otherwise load file.
    cache : :class:`bool`, (default True)
        Use cache of tiles data.
    """
    global _tiles

    if tilesfile is None:
        tilesfile = 'desi-tiles.fits'

    # - Check if tilesfile includes a path (absolute or relative)
    tilespath, filename = os.path.split(tilesfile)
    if tilespath == '':
        tilesfile = os.path.join(os.environ['LVMMODEL'], 'data', 'footprint', filename)

    # - standarize path location
    tilesfile = os.path.abspath(tilesfile)

    if cache and tilesfile in _tiles:
        tiledata = _tiles[tilesfile]
    else:
        with fits.open(tilesfile, memmap=False) as hdulist:
            tiledata = hdulist[1].data
        #
        # Temporary workaround for problem identified in
        # https://github.com/desihub/lvmmodel/issues/30
        #
        if any([c.bzero is not None for c in tiledata.columns]):
            foo = [_tiles[k].dtype for k in tiledata.dtype.names]

        # - Check for out-of-date tiles file
        if np.issubdtype(tiledata['OBSCONDITIONS'].dtype, 'u2'):
            import warnings
            warnings.warn('old desi-tiles.fits with uint16 OBSCONDITIONS; please update your $LVMMODEL checkout', DeprecationWarning)

        # - load cache for next time
        if cache:
            _tiles[tilesfile] = tiledata

    # - Filter to only the DESI footprint if requested
    subset = np.ones(len(tiledata), dtype=bool)
    if onlydesi:
        subset &= tiledata['IN_DESI'] > 0

    # - Filter out PROGRAM=EXTRA tiles if requested
    if not extra:
        subset &= ~np.char.startswith(tiledata['PROGRAM'], 'EXTRA')

    if np.all(subset):
        return tiledata
    else:
        return tiledata[subset]


_platescale = None


def load_platescale():
    '''
    Loads platescale.txt, returning structured array with columns

        radius: radius from center of focal plane [mm]
        theta: radial angle that has a centroid at this radius [deg]
        radial_platescale: Meridional (radial) plate scale [um/arcsec]
        az_platescale: Sagittal (azimuthal) plate scale [um/arcsec]
    '''
    global _platescale
    if _platescale is not None:
        return _platescale

    infile = findfile('focalplane/platescale.txt')
    columns = [
        ('radius', 'f8'),
        ('theta', 'f8'),
        ('radial_platescale', 'f8'),
        ('az_platescale', 'f8'),
    ]
    _platescale = np.loadtxt(infile, usecols=[0, 1, 6, 7], dtype=columns)
    return _platescale


def reset_cache():
    '''Reset I/O cache'''
    global _thru, _psf, _params, _gfa, _fiberpos, _tiles, _platescale
    _thru = dict()
    _psf = dict()
    _params = None
    _gfa = None
    _fiberpos = None
    _tiles = dict()
    _platescale = None


def load_target_info():
    '''
    Loads data/targets/targets.yaml and returns the nested dictionary

    This is primarily syntactic sugar to avoid end users constructing
    paths and filenames by hand (which e.g. broke when targets.dat was
    renamed to targets.yaml)
    '''
    targetsfile = os.path.join(datadir(), 'targets', 'targets.yaml')
    if not os.path.exists(targetsfile):
        targetsfile = os.path.join(datadir(), 'targets', 'targets.dat')

    with open(targetsfile) as fx:
        data = yaml.load(fx)

    return data


def load_pixweight(nside):
    '''
    Loads lvmmodel/data/footprint/desi-healpix-weights.fits

        nside: after loading, the array will be resampled to the
               passed HEALPix nside
    '''
    import healpy as hp
    # ADM read in the standard pixel weights file
    pixfile = os.path.join(os.environ['LVMMODEL'], 'data', 'footprint', 'desi-healpix-weights.fits')
    with fits.open(pixfile) as hdulist:
        pix = hdulist[0].data

    # ADM determine the file's nside, and flag a warning if the passed nside exceeds it
    npix = len(pix)
    truenside = hp.npix2nside(len(pix))
    if truenside < nside:
        log.warning("downsampling is fuzzy...Passed nside={}, "
                    "but file {} is stored at nside={}".format(nside, pixfile, truenside))

    # ADM resample the map
    return hp.pixelfunc.ud_grade(pix, nside, order_in='NESTED', order_out='NESTED')


def findfile(filename):
    '''
    Return full path to data file $LVMMODEL/data/filename

    Note: this is a precursor for a potential future refactor where
    lvmmodel data would be installed with the package and $LVMMODEL
    would become an optional override.
    '''
    return os.path.join(datadir(), filename)


def datadir():
    '''
    Returns location to lvmmodel data

    if set, $LVMMODEL overrides data installed with the package
    '''
    if 'LVMMODEL' in os.environ:
        return os.path.abspath(os.path.join(os.environ['LVMMODEL'], 'data'))
    else:
        import pkg_resources
        return pkg_resources.resource_filename('lvmmodel', 'data')
