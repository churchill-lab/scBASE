# -*- coding: utf-8 -*-

"""Main module."""

from __future__ import print_function
from . import utils
from . import get_data
from past.builtins import xrange
import os
import glob
import numpy as np
from scipy.sparse import dok_matrix
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
        num_genes, num_cells = ds.shape
        ds.ra['pi_P'] = dok_matrix((num_genes, 1), np.float64)
        ds.ra['pi_B'] = dok_matrix((num_genes, 1), np.float64)
        ds.ra['pi_M'] = dok_matrix((num_genes, 1), np.float64)
        ds.ra['alpha_mono'] = dok_matrix((num_genes, 1), np.float64)
        ds.ra['alpha_1'] = dok_matrix((num_genes, 1), np.float64)
        ds.ra['alpha_2'] = dok_matrix((num_genes, 1), np.float64)
        ds.ra['Rhat'] = dok_matrix((num_genes, 1), np.float64)
        ds.layers['pi_pk'] = 'float64'
        ds.layers['pi_bk'] = 'float64'
        ds.layers['pi_mk'] = 'float64'
        ds.layers['p_gk']  = 'float64'

        for f in flist:
            try:
                curdata_fh = np.load(f)['group1']
                curdata = curdata_fh.item()
            except IndexError:
                curdata_fh = np.load(f)
                curdata = curdata_fh.item()
            except IOError:
                LOG.warn("Error loading %s" % f)
            for g_key, g_fitting in curdata.iteritems():
                cur_gid = gid[g_key]
                ds.ra['pi_M'][cur_gid] = g_fitting['ase'][0, 0]
                ds.ra['pi_P'][cur_gid] = g_fitting['ase'][1, 0]
                ds.ra['pi_B'][cur_gid] = g_fitting['ase'][2, 0]
                ds.ra['alpha_1'][cur_gid] = g_fitting['ase'][4, 0]
                ds.ra['alpha_2'][cur_gid] = g_fitting['ase'][5, 0]
                ds.ra['Rhat'][cur_gid] = g_fitting['ase'][-1, -1]
                # Get ASE point estimation
                pi_gk = g_fitting['ase'][6+num_cells:6+num_cells*4, 0].reshape(3, num_cells)
                cur_theta = np.zeros(shape=pi_gk.shape)
                alpha_mono = g_fitting['ase'][3, 0]
                ds.ra['alpha_mono'][cur_gid] = alpha_mono
                cur_theta[0] = alpha_mono/(alpha_mono+1)
                cur_theta[1] = 1/(alpha_mono+1)
                cur_theta[2] = g_fitting['ase'][6:6+num_cells, 0]  # theta_{b,k}
                ds.layers['p_gk'][cur_gid, :] = (pi_gk * cur_theta).sum(axis=0)
        ds.close()
    elif filetype == "counts":
        pass



