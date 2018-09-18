# -*- coding: utf-8 -*-

"""Main module."""

from __future__ import print_function
from . import utils
from . import get_data
from past.builtins import xrange
import os
import glob
import numpy as np
import loompy
import pystan
try:
    import cPickle as pickle
except:
    import pickle

LOG = utils.get_logger()

try:
    xrange
except NameError:
    xrange = range


def disambiguate(alntools_file):
    LOG.warn('Quantifying allele-specific expression in each cell')
    LOG.info('Level-1 verbose')
    LOG.debug('Level-2 verbose')
    print(get_data('README'))


def __mcmc(x, n, model):
    data = {'N': len(n), 'n': n.astype('int'), 'x': x.astype('int')}
    fit_ase = model.sampling(data=data)
    LOG.debug(fit_ase)
    return fit_ase


def run_mcmc(loomfile, model, g_start, g_end):
    LOG.warn('Quantifying allele-specific expression in each cell')
    LOG.info('Level-1 verbose is on')
    LOG.debug('Level-2 verbose is also on')
    model_file = '%s.pkl' % model
    LOG.warn('Input file: %s' % get_data(model_file))
    stan_model = pickle.load(open(get_data(model_file), 'rb'))
    LOG.debug(print(stan_model.model_code))
    ds = loompy.connect(loomfile)
    if g_end < 0:
        g_end = ds.shape[0]
    elif g_end == 0:
        g_end = g_start
    outbase = 'scbase.%d-%d' % (g_start, g_end)
    param = dict()
    g_processed = 0
    tot_layer = ds.layers.keys()[0]
    mat_layer = ds.layers.keys()[1]
    for g in xrange(g_start, g_end):
        n = ds.layers[tot_layer][g]
        x = ds.layers[mat_layer][g]
        param[ds.row_attrs['gname'][g]] = __mcmc(x, n, stan_model)
        g_processed += 1
    LOG.info("All {:,d} genes have been processed.".format(g_processed))
    ds.close()
    np.savez_compressed('%s.%s' % (outbase, 'param.npz'), **param)


def collate(indir, loomfile, filetype, filename):
    if filename is not None:
        flist = glob.glob(os.path.join(indir, filename))
    else:
        if filetype == 'params':
            flist = glob.glob(os.path.join(indir, '*.param.npz'))
            LOG.info(os.path.join(indir, '*.param.npz'))
        elif filetype == 'counts':
            flist = glob.glob(os.path.join(indir, '*genes*counts'))
            LOG.info(os.path.join(indir, '*genes*counts'))
    LOG.info('Found %d files' % len(flist))
    if filetype == 'params':
        ds = loompy.connect(loomfile)
        gid = dict(zip(ds.row_attrs['gname'], np.arange(ds.shape[0])))


