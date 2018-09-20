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
    LOG.debug('ASE model code\n%s' % stan_model_ase.model_code)
    LOG.warn('TGX model file: %s' % get_data(model_file_tgx))
    stan_model_tgx = pickle.load(open(get_data(model_file_tgx), 'rb'))
    LOG.debug('TGX model code\n%s' % stan_model_tgx.model_code)
    ds = loompy.connect(loomfile)
    if end is None:
        end = start+1
    LOG.warn('Genes from %d to %d (0-based indexing)' % (start, end))
    libsz = ds.ca['size']
    c = libsz / np.median(libsz)
    LOG.debug('c: %s' % '\t'.join(c[:6].astype(str)))
    param = dict()
    processed = 0
    tgx_layer = ''
    mat_layer = hapcode[0]
    for g in xrange(start, end):
        if ds.ra['selected'][g]:
            LOG.warn('Loading data for Gene %s' % ds.ra['GeneID'][g])
            n = ds.layers[tgx_layer][g]
            x = ds.layers[mat_layer][g]
            LOG.debug('x: %s' % '\t'.join(x[:6].astype(int).astype(str)))
            LOG.debug('n: %s' % '\t'.join(n[:6].astype(int).astype(str)))
            cur_param = dict()
            LOG.warn('Fitting ASE with %s model' % model[0])
            cur_param['ase'] = __mcmc_ase(x, n, stan_model_ase).summary()['summary']
            LOG.warn('Fitting TGX with %s model' % model[1])
            cur_param['tgx'] = __mcmc_tgx(n, c, stan_model_tgx).summary()['summary']
            param[ds.row_attrs['GeneID'][g]] = cur_param
            processed += 1
    LOG.info("All {:,d} genes have been processed.".format(processed))
    outfile = os.path.join(outdir, 'scbase.%d-%d.param.npz' % (start, end))
    np.savez_compressed(outfile, **param)
    ds.close()


def collate(indir, loomfile, filetype, filename, model):
    if filename is not None:
        flist = glob.glob(os.path.join(indir, filename))
    else:
        if filetype == 'params':
            flist = glob.glob(os.path.join(indir, '*.param.npz'))
            LOG.info(os.path.join(indir, '*.param.npz'))
        elif filetype == 'counts':
            flist = glob.glob(os.path.join(indir, '*genes*counts'))
            LOG.info(os.path.join(indir, '*genes*counts'))
    LOG.info('Found %d file(s)' % len(flist))

    if filetype == 'params':
        ds = loompy.connect(loomfile)
        gid = dict(zip(ds.row_attrs['GeneID'], np.arange(ds.shape[0])))
        num_genes, num_cells = ds.shape

        # Initialize storage for ASE results
        if model[0] == 'zoibb':
            pi_p = dok_matrix((num_genes, 1), np.float64)
            pi_b = dok_matrix((num_genes, 1), np.float64)
            pi_m = dok_matrix((num_genes, 1), np.float64)
            alpha_mono = dok_matrix((num_genes, 1), np.float64)
            alpha_ase1 = dok_matrix((num_genes, 1), np.float64)
            alpha_ase2 = dok_matrix((num_genes, 1), np.float64)
            rhat_ase = dok_matrix((num_genes, 1), np.float64)
            ds.layers['pi_pk'] = 'float64'
            ds.layers['pi_bk'] = 'float64'
            ds.layers['pi_mk'] = 'float64'
            ds.layers['p_k'] = 'float64'
        else:
            raise NotImplementedError  # Add initiation for new ASE models here!!

        # Initialize storage for TGX results
        if model[1] == 'pg':
            alpha_tgx1 = dok_matrix((num_genes, 1), np.float64)
            alpha_tgx2 = dok_matrix((num_genes, 1), np.float64)
            rhat_tgx = dok_matrix((num_genes, 1), np.float64)
            ds.layers['lambda_k'] = 'float64'
        else:
            raise NotImplementedError  # Add initiation for new TGX models here!!

        for f in flist:
            LOG.warn('Loading %s' % f)
            curdata_fh = np.load(f)
            #curdata = curdata_fh.item()
            for g_key, g_results in curdata_fh.items():
                cur_gid = gid[g_key]
                LOG.debug('Current gene index: %d' % cur_gid)
                g_fitting = g_results.item()
                LOG.warn('Storing the fitting results of %s' % g_key)

                # Process ASE results
                if model[0] == 'zoibb':
                    LOG.info('Writing the ASE results by ZOIBB model')
                    pi_m[cur_gid] = g_fitting['ase'][0, 0]
                    pi_p[cur_gid] = g_fitting['ase'][1, 0]
                    pi_b[cur_gid] = g_fitting['ase'][2, 0]
                    LOG.debug('[ pi_p, pi_b, pi_m ] = [ %.3f %.3f %.3f ]' %
                              (g_fitting['ase'][1, 0], g_fitting['ase'][2, 0], g_fitting['ase'][0, 0]))
                    alpha_ase1[cur_gid] = g_fitting['ase'][4, 0]
                    alpha_ase2[cur_gid] = g_fitting['ase'][5, 0]
                    LOG.debug('[ alpha_ase1, alpha_ase2 ] = [ %.3f %.3f ]' %
                              (g_fitting['ase'][4, 0], g_fitting['ase'][5, 0]))
                    rhat_ase[cur_gid] = g_fitting['ase'][-1, -1]
                    LOG.debug('Rhat_ase = %.3f' % g_fitting['ase'][-1, -1])
                    # Get ASE point estimation
                    pi_k = g_fitting['ase'][6+num_cells:6+num_cells*4, 0].reshape(3, num_cells)
                    ds.layers['pi_mk'][cur_gid, :] = pi_k[0]
                    ds.layers['pi_pk'][cur_gid, :] = pi_k[1]
                    ds.layers['pi_bk'][cur_gid, :] = pi_k[2]
                    cur_theta = np.zeros(shape=pi_k.shape)
                    cur_alpha_mono = g_fitting['ase'][3, 0]
                    alpha_mono[cur_gid] = cur_alpha_mono
                    LOG.debug('alpha_mono = %.3f' % g_fitting['ase'][3, 0])
                    cur_theta[0] = cur_alpha_mono/(cur_alpha_mono+1)
                    cur_theta[1] = 1/(cur_alpha_mono+1)
                    cur_theta[2] = g_fitting['ase'][6:6+num_cells, 0]  # theta_{b,k}
                    ds.layers['p_k'][cur_gid, :] = (pi_k * cur_theta).sum(axis=0)
                else:  # Add handling of new ASE models here!!
                    raise NotImplementedError

                # Process TGX results
                if model[1] == 'pg':
                    LOG.info('Writing the TGX results by PG model')
                    # Get TGX point estimation
                    alpha_tgx1[cur_gid] = g_fitting['tgx'][0, 0]  # two alphas
                    alpha_tgx2[cur_gid] = g_fitting['tgx'][1, 0]  # two alphas
                    LOG.debug('[ alpha_tgx1, alpha_tgx2 ] = [ %.3f %.3f ]' %
                              (g_fitting['tgx'][4, 0], g_fitting['tgx'][5, 0]))
                    ds.layers['lambda_k'][cur_gid, :]  = g_fitting['tgx'][2:2+num_cells, 0]  # lambda_k
                    rhat_tgx[cur_gid] = g_fitting['tgx'][-1, -1]
                    LOG.debug('Rhat_tgx = %.3f' % g_fitting['tgx'][-1, -1])
                else:  # Add handling of new TGX models here!!
                    raise NotImplementedError

        # Store ASE results
        if model[0] == 'zoibb':
            ds.ra['pi_p'] = pi_p
            ds.ra['pi_b'] = pi_b
            ds.ra['pi_m'] = pi_m
            ds.ra['alpha_mono'] = alpha_mono
            ds.ra['alpha_ase1'] = alpha_ase1
            ds.ra['alpha_ase2'] = alpha_ase2
            ds.ra['Rhat_ase'] = rhat_ase
        else:  # Add storing of new ASE models here!!
            raise NotImplementedError

        # Store TGX results
        if model[1] == 'pg':
            ds.ra['alpha_tgx1'] = alpha_tgx1
            ds.ra['alpha_tgx2'] = alpha_tgx2
            ds.ra['Rhat_tgx'] = rhat_tgx
        else:  # Add storing of new TGX models here!!
            raise NotImplementedError
        ds.close()

    elif filetype == "counts":
        pass



