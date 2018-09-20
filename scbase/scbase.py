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


def disambiguate(alntools_file, start, end):
    LOG.warn('Quantifying allele-specific expression in each cell')
    LOG.info('Level-1 verbose')
    LOG.debug('Level-2 verbose')
    print(get_data('README'))


def __mcmc_ase(x, n, model):
    data = {'N': len(n), 'n': n.astype('int'), 'x': x.astype('int')}
    fit_ase = model.sampling(data=data)
    LOG.debug(fit_ase)
    return fit_ase


def __mcmc_tgx(n, c, model):
    data = {'N': len(n), 'n': n.astype('int'), 'C':c}
    fit_tgx = model.sampling(data=data)
    LOG.debug(fit_tgx)
    return fit_tgx


def run_mcmc(loomfile, model, hapcode, start, end, outdir):
    LOG.warn('Quantifying allele-specific expression in each cell')
    LOG.info('Level-1 verbose is on')
    LOG.debug('Level-2 verbose is also on')
    model_file_ase = '%s.pkl' % model[0]
    model_file_tgx = '%s.pkl' % model[1]
    LOG.warn('ASE model file: %s' % get_data(model_file_ase))
    stan_model_ase = pickle.load(open(get_data(model_file_ase), 'rb'))
    LOG.debug(stan_model_ase.model_code)
    LOG.warn('TGX model file: %s' % get_data(model_file_tgx))
    stan_model_tgx = pickle.load(open(get_data(model_file_tgx), 'rb'))
    LOG.debug(stan_model_tgx.model_code)
    ds = loompy.connect(loomfile)
    if end is None:
        end = start+1
    LOG.warn('Genes from %d to %d (0-based indexing)' % (start, end))
    libsz = ds.ca['size']
    c = libsz / np.median(libsz)
    param = dict()
    processed = 0
    tgx_layer = ''
    mat_layer = hapcode[0]
    for g in xrange(start, end):
        if ds.ra['selected'][g]:
            LOG.warn('Loading data for Gene %s [%s]' % (ds.ra['gname'][g], ds.ra['gsymb'][g]))
            n = ds.layers[tgx_layer][g]
            x = ds.layers[mat_layer][g]
            LOG.debug()
            cur_param = dict()
            LOG.warn('Fitting ASE with %s model' % model[0])
            cur_param['ase'] = __mcmc_ase(x, n, stan_model_ase)
            LOG.warn('Fitting TGX with %s model' % model[1])
            cur_param['tot'] = __mcmc_tgx(n, c, stan_model_tgx)
            param[ds.row_attrs['gname'][g]] = cur_param
            processed += 1
    LOG.info("All {:,d} genes have been processed.".format(processed))
    outfile = os.path.join(outdir, 'scbase.%d-%d' % (start, end), 'param.npz')
    np.savez_compressed(outfile, **param)
    ds.close()


def collate(indir, loomfile, filetype, filename, groupname, model):
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
        if model[0] == 'zoibb':
            ds.ra['pi_P'] = dok_matrix((num_genes, 1), np.float64)
            ds.ra['pi_B'] = dok_matrix((num_genes, 1), np.float64)
            ds.ra['pi_M'] = dok_matrix((num_genes, 1), np.float64)
            ds.ra['alpha_mono'] = dok_matrix((num_genes, 1), np.float64)
            ds.ra['alpha_ase1'] = dok_matrix((num_genes, 1), np.float64)
            ds.ra['alpha_ase2'] = dok_matrix((num_genes, 1), np.float64)
            ds.ra['Rhat_ase'] = dok_matrix((num_genes, 1), np.float64)
            ds.layers['pi_pk'] = 'float64'
            ds.layers['pi_bk'] = 'float64'
            ds.layers['pi_mk'] = 'float64'
            ds.layers['p_gk']  = 'float64'
        else:
            raise NotImplementedError  # Add initiation for new ASE models here!!
        if model[1] == 'pg':
            ds.ra['alpha_tgx1'] = dok_matrix((num_genes, 1), np.float64)
            ds.ra['alpha_tgx2'] = dok_matrix((num_genes, 1), np.float64)
            ds.ra['Rhat_tgx'] = dok_matrix((num_genes, 1), np.float64)
            ds.layers['lambda_gk']  = 'float64'
        else:
            raise NotImplementedError  # Add initiation for new TGX models here!!

        for f in flist:
            try:
                curdata_fh = np.load(f)[groupname]
                curdata = curdata_fh.item()
            except IndexError:
                curdata_fh = np.load(f)
                curdata = curdata_fh.item()
            except IOError:
                LOG.warn("Error loading %s" % f)
            for g_key, g_fitting in curdata.iteritems():
                cur_gid = gid[g_key]
                if model[0] == 'zoibb':
                    ds.ra['pi_M'][cur_gid] = g_fitting['ase'][0, 0]
                    ds.ra['pi_P'][cur_gid] = g_fitting['ase'][1, 0]
                    ds.ra['pi_B'][cur_gid] = g_fitting['ase'][2, 0]
                    ds.ra['alpha_ase1'][cur_gid] = g_fitting['ase'][4, 0]
                    ds.ra['alpha_ase2'][cur_gid] = g_fitting['ase'][5, 0]
                    ds.ra['Rhat_ase'][cur_gid] = g_fitting['ase'][-1, -1]
                    # Get ASE point estimation
                    pi_gk = g_fitting['ase'][6+num_cells:6+num_cells*4, 0].reshape(3, num_cells)
                    cur_theta = np.zeros(shape=pi_gk.shape)
                    alpha_mono = g_fitting['ase'][3, 0]
                    ds.ra['alpha_mono'][cur_gid] = alpha_mono
                    cur_theta[0] = alpha_mono/(alpha_mono+1)
                    cur_theta[1] = 1/(alpha_mono+1)
                    cur_theta[2] = g_fitting['ase'][6:6+num_cells, 0]  # theta_{b,k}
                    ds.layers['p_gk'][cur_gid, :] = (pi_gk * cur_theta).sum(axis=0)
                else:
                    raise NotImplementedError  # Add handling of new ASE models here!!
                # Get TGX point estimation
                if model[1] == 'pg':
                    ds.ra['alpha_tgx1'][cur_gid] = g_fitting['tot'][0, 0]  # two alphas
                    ds.ra['alpha_tgx2'][cur_gid] = g_fitting['tot'][1, 0]  # two alphas
                    ds.layers['lambda_gk'][cur_gid, :]  = g_fitting['tot'][2:2+num_cells, 0]  # lambda_gk
                    ds.ra['Rhat_tgx'][cur_gid] = g_fitting['tot'][-1, -1]
                else:
                    raise NotImplementedError  # Add handling of new TGX models here!!
        ds.close()

    elif filetype == "counts":
        pass



