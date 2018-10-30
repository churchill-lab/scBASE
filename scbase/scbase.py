# -*- coding: utf-8 -*-

"""Main module."""

from __future__ import print_function
from . import utils
from . import get_data
from past.builtins import xrange
import os
import time
import glob
import numpy as np
from scipy.sparse import lil_matrix, dok_matrix, csr_matrix, hstack
from scipy.special import digamma, polygamma
import loompy
import pystan
from subprocess import call
try:
    import cPickle as pickle
except:
    import pickle

LOG = utils.get_logger()

try:
    xrange
except NameError:
    xrange = range


# def disambiguate(alnfile, start, end):
#     LOG.warn('Quantifying allele-specific expression in each cell')
#     LOG.info('Verbose level 1 [ON]')
#     LOG.debug('Verbose level 2 [ON]')
#     raise NotImplementedError('Coming soon once alntools and emase-zero projects are completed.')


def select(loomfile, min_read_count, min_cell_count, layer):
    with loompy.connect(loomfile) as ds:
        gsurv = (ds.sparse(layer=layer) >= min_read_count).sum(axis=1) > min_cell_count
        ds.ra.Selected = np.squeeze(np.asarray(gsurv))
        LOG.info('Total %d genes selected' % gsurv.sum())
        # totals = ds.map([np.sum], axis=1)[0]  # Select based upon cell size?


def __mcmc4ase(x, n, model):
    data = {'N': len(n), 'n': n.astype('int'), 'x': x.astype('int')}
    fit_ase = model.sampling(data=data)
    LOG.debug(fit_ase)
    return fit_ase


def __mcmc4tgx(n, c, model):
    data = {'N': len(n), 'n': n.astype('int'), 'C':c}
    fit_tgx = model.sampling(data=data)
    LOG.debug(fit_tgx)
    return fit_tgx


def run_mcmc(loomfile, model, hapcode, start, end, outfile):
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
    ds = loompy.connect(loomfile, 'r')
    if end is None:
        end = ds.shape[0]
    LOG.warn('Genes from %d to %d (0-based indexing)' % (start, end))
    c = ds.ca.Size / np.median(ds.ca.Size)
    LOG.debug('c: %s' % '\t'.join(c[:6].astype(str)))
    param = dict()
    processed = 0
    tgx_layer = ''
    mat_layer = hapcode[0]
    for g in xrange(start, end):
        if ds.ra.Selected[g]:
            LOG.warn('Loading data for Gene %s' % ds.ra['GeneID'][g])
            n = ds.layers[tgx_layer][g]
            x = ds.layers[mat_layer][g]
            LOG.debug('x: %s' % '\t'.join(x[:6].astype(int).astype(str)))
            LOG.debug('n: %s' % '\t'.join(n[:6].astype(int).astype(str)))
            cur_param = dict()
            LOG.warn('Fitting ASE with %s model' % model[0])
            cur_param['ase'] = __mcmc4ase(x, n, stan_model_ase).summary()['summary']
            LOG.warn('Fitting TGX with %s model' % model[1])
            cur_param['tgx'] = __mcmc4tgx(n, c, stan_model_tgx).summary()['summary']
            param[ds.row_attrs['GeneID'][g]] = cur_param
            processed += 1
    LOG.info("All {:,d} genes have been processed.".format(processed))
    if outfile is None:
        outfile = 'scbase.%05d-%05d.param.npz' % (start, end)
    np.savez_compressed(outfile, **param)
    ds.close()


def run_mcmc_from_npz(datafile, model, hapcode, start, end, outfile):
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
    data_dict = np.load(datafile)
    if end is None:
        end = data_dict['shape'][0]
    LOG.warn('Genes from %d to %d (0-based indexing)' % (start, end))
    libsz = data_dict['Size']
    c = libsz / np.median(libsz)
    LOG.debug('c: %s' % '\t'.join(c[:6].astype(str)))
    param = dict()
    processed = 0
    tgx_layer = ''
    mat_layer = hapcode[0]
    dmat_dict = data_dict['Counts'].item()
    for g in xrange(start, end):
        if data_dict['Selected'][g]:
            LOG.warn('Loading data for Gene %s' % data_dict['GeneID'][g])
            n = dmat_dict[tgx_layer][g]
            x = dmat_dict[mat_layer][g]
            LOG.debug('x: %s' % '\t'.join(x[:6].astype(int).astype(str)))
            LOG.debug('n: %s' % '\t'.join(n[:6].astype(int).astype(str)))
            cur_param = dict()
            LOG.warn('Fitting ASE with %s model' % model[0])
            cur_param['ase'] = __mcmc4ase(x, n, stan_model_ase).summary()['summary']
            LOG.warn('Fitting TGX with %s model' % model[1])
            cur_param['tgx'] = __mcmc4tgx(n, c, stan_model_tgx).summary()['summary']
            param[data_dict['GeneID'][g]] = cur_param
            processed += 1
    LOG.info("All {:,d} genes have been processed.".format(processed))
    if outfile is None:
        outfile = 'scbase.%05d-%05d.param.npz' % (start, end)
    np.savez_compressed(outfile, **param)


def __update_shape(sufficient, tol=0.000001, max_iters=100):
    shape = 0.5 / sufficient
    for cur_iter in range(max_iters):
        shape_prev = shape.copy()
        g = np.log(shape) - sufficient - digamma(shape)
        h = 1/shape - polygamma(1, shape)
        shape = 1 / (1/shape + g/(np.power(shape, 2) * h))
        abs_err = np.abs(shape - shape_prev)
        abs_err_max = np.max(abs_err)
        abs_err_arg = np.argmax(abs_err)
        if abs_err_max < tol:
            LOG.debug("__update_shape ran total of %d iterations (Max error=%.9f)" % (cur_iter+1, abs_err_max))
            break
        LOG.debug('Iter #%04d: Error = %.4f [argmax = %d]' % (cur_iter, abs_err_max, abs_err_arg))
    return shape


def __em4tgx(cntmat, scaler, percentile, tol, max_iters):
    cntmat_scaled = cntmat.copy()
    cntmat_scaled.data = cntmat_scaled.data / scaler[cntmat_scaled.indices]
    # libsz_scaled = np.squeeze(np.asarray(cntmat_scaled.sum(axis=0)))
    ridx = np.repeat(np.arange(cntmat.shape[0]), np.diff(cntmat.indptr))
    #
    # Initial mu and phi
    #
    x_nnz = np.squeeze(np.asarray((cntmat > 0).sum(axis=1)))
    mean_x_scaled = np.squeeze(np.asarray(cntmat_scaled.sum(axis=1))) / x_nnz
    mean_x_scaled_square = np.squeeze(np.asarray(cntmat_scaled.power(2).sum(axis=1))) / x_nnz
    var_x_scaled = mean_x_scaled_square - np.power(mean_x_scaled, 2)
    phi = var_x_scaled / mean_x_scaled - 1
    phi[phi < 0] = 0.001
    mu = mean_x_scaled / phi
    p = cntmat.copy()
    lamb = cntmat.copy()
    log_lamb = cntmat.copy()

    for cur_iter in range(max_iters):
        #
        # Initialize an iteration
        #
        mu_prev = mu.copy()

        #
        # E-step
        #
        p.data = phi[ridx] / (scaler[cntmat.indices]*phi[ridx] + 1)
        x_plus_mu = cntmat.copy()
        x_plus_mu.data += mu[ridx]
        lamb.data = x_plus_mu.data * p.data
        mean_lamb = np.squeeze(np.asarray(lamb.sum(axis=1))) / x_nnz
        log_lamb.data = digamma(x_plus_mu.data) + np.log(p.data)
        mean_log_lamb = np.squeeze(np.asarray(log_lamb.sum(axis=1))) / x_nnz
        suff = np.log(mean_lamb) - mean_log_lamb

        #
        # M-step
        #
        mu = __update_shape(suff, tol, max_iters)
        phi = mean_lamb / mu

        #
        # Check termination
        #
        err = np.abs(mu - mu_prev)
        err_pct = np.percentile(err, percentile)
        num_converged = sum(err < tol)
        if err_pct < tol:
            break
        LOG.warn('Iter #%04d: %6s genes converged below the tolerance level of %.1E' % (cur_iter+1, num_converged, tol))
        LOG.debug('Median error=%.6f' % err_pct)
    if cur_iter+1 == max_iters:
        LOG.warn('Reached the maximum number of iterations')
    return lamb, mu, phi, err


def __em4ase(cntmat, tol, max_iters):
    raise NotImplementedError('EM algorithm for ASE is coming soon.')


def run_em(loomfile, model, common_scale, percentile, hapcode, start, end, tol, max_iters, outfile):
    if model[0] == 'null' and model[1] == 'null':
        raise RuntimeError('At least either of ASE or TGX model should be specified.')
    # ASE model
    if model[0] == 'null':
        LOG.warn('No ASE model will run.')
    elif model[0] == 'zoibb':
        raise NotImplementedError('EM version of ZOIBB model is coming soon')
    else:
        raise NotImplementedError('Only ZOIBB model will be available for ASE in run_em')
    # TGX model
    if model[1] == 'null':
        LOG.warn('No TGX model will run.')
    elif model[1] == 'pg':
        with loompy.connect(loomfile) as ds:
            num_genes, num_cells = ds.shape
            LOG.info('Loading data from %s' % loomfile)
            origmat = ds.sparse().tocsr()
            LOG.info('Processing data matrix')
            if 'Selected' in ds.ca.keys():
                csurv = np.where(ds.ca.Selected > 0)[0]
                cntmat = origmat[:, csurv]
            else:
                csurv = np.arange(num_cells)
                cntmat = origmat
            LOG.info('The number of selected cells: %d' % len(csurv))
            
            libsz = np.squeeze(np.asarray(cntmat.sum(axis=0)))
            scaler = libsz / common_scale
            if 'Selected' in ds.ra.keys():
                gsurv1 = ds.ra.Selected > 0
            else:
                gsurv1 = np.ones(num_genes)
            gsurv2 = np.squeeze(np.asarray((cntmat > 0).sum(axis=1) > 0))
            gsurv = np.where(np.logical_and(gsurv1, gsurv2))[0]
            LOG.info('The number of selected genes: %d' % len(gsurv))
            cntmat = cntmat[gsurv, :]
            LOG.info('Running EM algorithm for TGX')
            lambda_mat, mu, phi, err = __em4tgx(cntmat, scaler, percentile, tol, max_iters)
            LOG.info('There were %d genes that converged below the tolerance level of %.1E' % (sum(err < tol), tol))
            LOG.info('Saving results to %s' % loomfile)
            resmat = csr_matrix((origmat.shape))
            resmat.indptr = np.ones(resmat.indptr.shape, dtype='int') * lambda_mat.indptr[-1]
            resmat.indptr[0] = 0
            resmat_indptr_part = np.repeat(lambda_mat.indptr[1:-1], np.diff(gsurv))
            resmat.indptr[1:len(resmat_indptr_part)+1] = resmat_indptr_part
            resmat.indices = csurv[lambda_mat.indices]
            resmat.data = lambda_mat.data
            ds.layers['lambda'] = resmat
            mu_res = dok_matrix((num_genes, 1), float)
            mu_res[gsurv] = mu[:, np.newaxis]
            ds.ra['mu'] = mu_res
            phi_res = dok_matrix((num_genes, 1), float)
            phi_res[gsurv] = phi[:, np.newaxis]
            ds.ra['phi'] = phi_res
            err_res = dok_matrix((num_genes, 1), float)
            err_res[gsurv] = err[:, np.newaxis]
            ds.ra['err'] = err_res
            g_selected = dok_matrix((num_genes, 1), float)
            g_selected[gsurv] = 1
            ds.ra['Selected:TGX:EM'] = g_selected
        LOG.info("Finished EM for TGX")
    else:
        raise NotImplementedError('Only Gamma-Poisson model is available for TGX in run_em.')


def submit(loomfile, model, hapcode, chunk, outdir, email, queue, mem, walltime, systype, dryrun):
    LOG.warn('Loom file: %s' % loomfile)
    LOG.warn('Models: %s, %s' % (model[0], model[1]))
    LOG.warn('HPC system type: %s' % systype)
    if dryrun:
        LOG.warn('Showing submission script only')

    with loompy.connect(loomfile, 'r') as ds:
        gsurv = np.where(ds.ra.Selected)[0]
        num_gsurv = len(gsurv)
        num_genes, num_cells = ds.shape
        LOG.warn('The number of selected genes: %d' % num_gsurv)
        LOG.warn('The number of selected cells: %d' % num_cells)
        LOG.warn('%d jobs will be submitted' % int(np.ceil(num_gsurv/chunk)))
    processed = 0

    if systype == 'pbs':
        tot_layer = ''
        mat_layer = hapcode[0]
        for idx_start in xrange(0, num_gsurv, chunk):
            idx_end = min(idx_start+chunk, num_gsurv-1)
            start = gsurv[idx_start]
            if idx_end < num_gsurv-1:
                end = gsurv[idx_end]
                genes = gsurv[idx_start:idx_end]
            else:  #idx_end == num_gsurv-1:
                end = num_genes
                genes = gsurv[idx_start:]
            LOG.info('Chunk start: %d, end %d' % (start, end))
            infile = os.path.join(outdir, '_chunk.%05d-%05d.npz' % (start, end))
            LOG.debug('Genes: %s' % ' '.join(genes.astype(str)))
            LOG.debug('Total %d genes submitted in this job' % len(genes))
            data_dict = dict()
            data_dict['shape'] = (len(genes), num_cells)
            with loompy.connect(loomfile, 'r') as ds:
                data_dict['GeneID'] = ds.ra.GeneID[genes]
                cur_chunk = dict()
                cur_chunk[tot_layer] = ds.layers[tot_layer][genes, :]
                cur_chunk[mat_layer] = ds.layers[mat_layer][genes, :]
                data_dict['Counts'] = cur_chunk
                data_dict['Size'] = ds.ca.Size
                data_dict['Selected'] = np.ones(len(genes))  # select all
                np.savez_compressed(infile, **data_dict)
            outfile = os.path.join(outdir, 'scbase.%05d-%05d.param.npz' % (start, end))
            job_par = 'ASE_MODEL=%s,TGX_MODEL=%s,MAT_HAPCODE=%s,PAT_HAPCODE=%s,OUTFILE=%s,INFILE=%s' % \
                      (model[0], model[1], hapcode[0], hapcode[1], outfile, infile)
            cmd = ['qsub']
            if email is not None:
                cmd += ['-M', email]
            if queue is not None:
                cmd += ['-q', queue]
            if mem > 0:
                cmd += ['-l', 'mem=%d' % mem]
            if walltime > 0:
                cmd += ['-l', 'walltime=%d:00:00' % walltime]
            cmd += ['-v', job_par]
            cmd += [os.path.join(os.path.dirname(os.environ['_']), 'run_mcmc_on_cluster.sh')]
            if dryrun:
                print(" ".join(cmd))
            else:
                LOG.info(" ".join(cmd))
                call(cmd)
                time.sleep(1.0)
            processed += len(genes)
        LOG.debug('Total %d genes were submitted' % processed)
        LOG.warn('Job submission complete')
    elif systype == 'pbs-with-whole-loom':  # Do not use this: loom is not stable
        for idx_start in xrange(0, num_gsurv, chunk):
            idx_end = min(idx_start+chunk, num_gsurv-1)
            start = gsurv[idx_start]
            if idx_end < num_gsurv-1:
                end = gsurv[idx_end]
                genes = gsurv[idx_start:idx_end]
            else:  #idx_end == num_gsurv-1:
                end = num_genes
                genes = gsurv[idx_start:]
            LOG.info('Chunk start: %d, end %d' % (start, end))
            LOG.debug('Genes: %s' % ' '.join(genes.astype(str)))
            LOG.debug('Total %d genes submitted in this job' % len(genes))
            outfile = os.path.join(outdir, 'scbase.%05d-%05d.param.npz' % (start, end))
            job_par = 'ASE_MODEL=%s,TGX_MODEL=%s,MAT_HAPCODE=%s,PAT_HAPCODE=%s,START=%d,END=%d,OUTFILE=%s,INFILE=%s' % \
                      (model[0], model[1], hapcode[0], hapcode[1], start, end, outfile, loomfile)
            cmd = ['qsub']
            if email is not None:
                cmd += ['-M', email]
            if queue is not None:
                cmd += ['-q', queue]
            if mem > 0:
                cmd += ['-l', 'mem=%d' % mem]
            if walltime > 0:
                cmd += ['-l', 'walltime=%d:00:00' % walltime]
            cmd += ['-v', job_par]
            cmd += [os.path.join(os.path.dirname(os.environ['_']), 'run_mcmc_on_cluster.sh')]
            if dryrun:
                print(" ".join(cmd))
            else:
                LOG.info(" ".join(cmd))
                call(cmd)
                time.sleep(1.0)
            processed += len(genes)
        LOG.debug('Total %d genes were submitted' % processed)
        LOG.warn('Job submission complete')
    elif systype == 'pbs-with-loom-chunks':  # Do not use this: loompy does not support this
        for idx_start in xrange(0, num_gsurv, chunk):
            idx_end = min(idx_start+chunk, num_gsurv-1)
            start = gsurv[idx_start]
            end = gsurv[idx_end]
            if idx_end < num_gsurv-1:
                end = gsurv[idx_end]
                genes = gsurv[idx_start:idx_end]
            else:  #idx_end == num_gsurv-1:
                end = num_genes
                genes = gsurv[idx_start:]
            LOG.info('Chunk start: %d, end %d' % (start, end))
            infile = os.path.join(outdir, '_chunk.%05d-%05d.loom' % (start, end))
            LOG.debug('Genes: %s' % ' '.join(genes.astype(str)))
            LOG.debug('Total %d genes submitted in this job' % len(genes))
            with loompy.connect(loomfile, 'r') as ds:
                with loompy.new(infile) as dsout:
                    for (ix, selection, view) in ds.scan(items=genes, axis=0):
                        LOG.debug('Genes in this view: %s' % ' '.join(selection.astype()))
                        dsout.add_columns(view.layers, col_attrs=view.col_attrs, row_attrs=view.row_attrs)
            outfile = os.path.join(outdir, 'scbase.%05d-%05d.param.npz' % (start, end))
            job_par = 'ASE_MODEL=%s,TGX_MODEL=%s,MAT_HAPCODE=%s,PAT_HAPCODE=%s,OUTFILE=%s,INFILE=%s' % \
                      (model[0], model[1], hapcode[0], hapcode[1], outfile, infile)
            cmd = ['qsub']
            if email is not None:
                cmd += ['-M', email]
            if queue is not None:
                cmd += ['-q', queue]
            if mem > 0:
                cmd += ['-l', 'mem=%d' % mem]
            if walltime > 0:
                cmd += ['-l', 'walltime=%d:00:00' % walltime]
            cmd += ['-v', job_par]
            cmd += [os.path.join(os.path.dirname(os.environ['_']), 'run_mcmc_on_cluster.sh')]
            if dryrun:
                print(" ".join(cmd))
            else:
                LOG.info(" ".join(cmd))
                call(cmd)
                time.sleep(1.0)
            processed += len(genes)
        LOG.debug('Total %d genes were submitted' % processed)
        LOG.warn('Job submission complete')
    elif 'lsf':
        raise NotImplementedError('LSF submission is not yet supported')
    else:
        raise RuntimeError('No plan to support other job scheduling system until we see many requests')


def collate(indir, loomfile, tidfile, filetype, filename, model):
    if model[0] == 'null' and model[1] == 'null':
        raise RuntimeError('At least either of ASE or TGX model should be specified.')

    if filetype == "counts":
        LOG.warn('Looking at %s directly for count files...' % os.path.abspath(indir))
        if filename is None:
            flist = glob.glob(os.path.join(indir, '*gene*counts'))
        else:
            flist = glob.glob(os.path.join(indir, filename))
        if len(flist) > 0:
            LOG.warn('%d files were found under %s' % (len(flist), indir))
            flist.sort()
        else:
            raise FileNotFoundError('No files to collate')

        f = flist[0]
        with open(f) as fh:
            curline = fh.readline()
            item = curline.rstrip().split('\t')
            hapcodes = item[1:-1]
            num_haps = len(hapcodes)
        LOG.warn('Haplotypes: %s' % '\t'.join(hapcodes))

        if tidfile is not None:
            geneID = np.loadtxt(tidfile, dtype=str, usecols=0)
            num_genes = len(geneID)
            gene_idx = dict(zip(geneID, np.arange(num_genes)))
        LOG.warn('Number of genes: %d [%s %s ...]' % (len(geneID), geneID[0], geneID[1]))

        dmat = dict()
        dmat[''] = lil_matrix((num_genes, 0))
        for h in hapcodes:
            dmat[h] = lil_matrix((num_genes, 0))
        
        cellID = list()
        # cix = -1
        for f in flist:
            new_data = np.zeros(0)
            with open(f) as fh:
                LOG.warn("Loading counts from %s:" % f)
                fh.readline()  # skip the header in each file
                for curline in fh:
                    item = curline.rstrip().split()
                    if '#sample_id' in curline:
                        if new_data.sum() > 0:
                            LOG.info("Storing results of Cell: %s" % cellID[-1])
                            dmat[''] = hstack((dmat[h], new_data[:, -1]))
                            for hix, h in enumerate(hapcodes):
                                dmat[h] = hstack((dmat[h], new_data[:, hix]))
                        new_data = lil_matrix((num_genes, num_haps+1))
                        cellID.append(item[1])
                        # cix += 1
                    else:
                        gi = gene_idx[item[0]]
                        if float(item[-1]) > 0:
                            new_data[gi] = np.array(item[1:]).astype(float)
            LOG.info("Storing results of Cell: %s" % cellID[-1])
            dmat[''] = hstack((dmat[h], new_data[:, -1]))
            for hix, h in enumerate(hapcodes):
                dmat[h] = hstack((dmat[h], new_data[:, hix]))
            LOG.info('All counts loaded from %s' % f)
        loompy.create(loomfile, dmat[''], row_attrs={'GeneID': geneID}, col_attrs={'CellID': np.array(cellID).astype(str)})
        LOG.warn('Created %s' % loomfile)
        ds = loompy.connect(loomfile)
        for h in hapcodes:
            ds.layers[h] = dmat[h]
        ds.ca['Size'] = dmat[''].sum(axis=0)        
        ds.close()
        LOG.warn('Done. You can add more row_attrs or col_attrs to %s' % loomfile)

    elif filetype == 'params':
        LOG.warn('Looking at %s directly for param files...' % os.path.abspath(indir))
        if filename is None:
            flist = glob.glob(os.path.join(indir, '*.param.npz'))
        else:
            flist = glob.glob(os.path.join(indir, filename))
        if len(flist) < 1:
            raise FileNotFoundError('No param files to collate')
        else:
            LOG.info('Found %d param file(s)' % len(flist))

        ds = loompy.connect(loomfile)
        gid = dict(zip(ds.row_attrs['GeneID'], np.arange(ds.shape[0])))
        num_genes, num_cells = ds.shape

        # Initialize storage for ASE results
        if model[0] == 'null':
            LOG.warn('ASE model will not run')
        elif model[0] == 'zoibb':
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
        # Add initiation for new ASE models here!!
        else:
            raise NotImplementedError('%s model does not exist' % model[0])

        # Initialize storage for TGX results
        if model[1] == 'null':
            LOG.warn('TGX model will not run')
        elif model[1] == 'pg':
            alpha_tgx1 = dok_matrix((num_genes, 1), np.float64)
            alpha_tgx2 = dok_matrix((num_genes, 1), np.float64)
            rhat_tgx = dok_matrix((num_genes, 1), np.float64)
            ds.layers['lambda_k'] = 'float64'
        # Add initiation for new TGX models here!!
        else:
            raise NotImplementedError('%s model does not exist' % model[1])

        for f in flist:
            LOG.info('Loading %s' % f)
            curdata_fh = np.load(f)
            for g_key, g_results in curdata_fh.items():
                cur_gid = gid[g_key]
                LOG.debug('Current gene index: %d' % cur_gid)
                g_fitting = g_results.item()
                LOG.warn('Storing the fitting results of %s' % g_key)

                # Process ASE results
                if model[0] == 'null':
                    pass
                elif model[0] == 'zoibb':
                    LOG.info('ASE results by ZOIBB model')
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
                # Add handling of new ASE models here!!
                else:
                    raise NotImplementedError('scBASE does not know how to process %s model results' % model[0])

                # Process TGX results
                if model[1] == 'null':
                    pass
                elif model[1] == 'pg':
                    LOG.info('TGX results by PG model')
                    # Get TGX point estimation
                    alpha_tgx1[cur_gid] = g_fitting['tgx'][0, 0]  # two alphas
                    alpha_tgx2[cur_gid] = g_fitting['tgx'][1, 0]  # two alphas
                    LOG.debug('[ alpha_tgx1, alpha_tgx2 ] = [ %.3f %.3f ]' %
                              (g_fitting['tgx'][4, 0], g_fitting['tgx'][5, 0]))
                    ds.layers['lambda_k'][cur_gid, :]  = g_fitting['tgx'][2:2+num_cells, 0]  # lambda_k
                    rhat_tgx[cur_gid] = g_fitting['tgx'][-1, -1]
                    LOG.debug('Rhat_tgx = %.3f' % g_fitting['tgx'][-1, -1])
                # Add handling of new TGX models here!!
                else:
                    raise NotImplementedError('scBASE does not know how to process %s model results' % model[1])
            LOG.warn('Finished processing %s' % f)

        # Store ASE results
        if model[0] == 'null':
            pass
        elif model[0] == 'zoibb':
            ds.ra['pi_p'] = pi_p
            ds.ra['pi_b'] = pi_b
            ds.ra['pi_m'] = pi_m
            ds.ra['alpha_mono'] = alpha_mono
            ds.ra['alpha_ase1'] = alpha_ase1
            ds.ra['alpha_ase2'] = alpha_ase2
            ds.ra['Rhat_ase'] = rhat_ase
        # Add storing of new ASE models here!!
        else:
            raise NotImplementedError('scBASE does not know how to store %s model results' % model[0])

        # Store TGX results
        if model[1] == 'null':
            pass
        elif model[1] == 'pg':
            ds.ra['alpha_tgx1'] = alpha_tgx1
            ds.ra['alpha_tgx2'] = alpha_tgx2
            ds.ra['Rhat_tgx'] = rhat_tgx
        # Add storing of new TGX models here!!
        else:
            raise NotImplementedError('scBASE does not know how to store %s model results' % model[1])
        ds.close()

    else:
        raise RuntimeError('filetype option should be either of --counts or --params')
