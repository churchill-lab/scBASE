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
from scipy.sparse import dok_matrix
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


def disambiguate(alnfile, start, end):
    LOG.warn('Quantifying allele-specific expression in each cell')
    LOG.info('Verbose level 1 [ON]')
    LOG.debug('Verbose level 2 [ON]')
    raise NotImplementedError('Coming soon once alntools and emase-zero projects are completed.')


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


def run_mcmc_from_loom(loomfile, model, hapcode, start, end, outfile):
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
        end = ds.shape[0]
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
    if outfile is None:
        outfile = 'scbase.%05d-%05d.param.npz' % (start, end)
    np.savez_compressed(outfile, **param)
    ds.close()


def run_mcmc(datafile, model, hapcode, start, end, outfile):
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
    libsz = data_dict['size']
    c = libsz / np.median(libsz)
    LOG.debug('c: %s' % '\t'.join(c[:6].astype(str)))
    param = dict()
    processed = 0
    tgx_layer = ''
    mat_layer = hapcode[0]
    dmat_dict = data_dict['counts'].item()
    for g in xrange(start, end):
        if data_dict['selected'][g]:
            LOG.warn('Loading data for Gene %s' % data_dict['GeneID'][g])
            n = dmat_dict[tgx_layer][g]
            x = dmat_dict[mat_layer][g]
            LOG.debug('x: %s' % '\t'.join(x[:6].astype(int).astype(str)))
            LOG.debug('n: %s' % '\t'.join(n[:6].astype(int).astype(str)))
            cur_param = dict()
            LOG.warn('Fitting ASE with %s model' % model[0])
            cur_param['ase'] = __mcmc_ase(x, n, stan_model_ase).summary()['summary']
            LOG.warn('Fitting TGX with %s model' % model[1])
            cur_param['tgx'] = __mcmc_tgx(n, c, stan_model_tgx).summary()['summary']
            param[data_dict['GeneID'][g]] = cur_param
            processed += 1
    LOG.info("All {:,d} genes have been processed.".format(processed))
    if outfile is None:
        outfile = 'scbase.%05d-%05d.param.npz' % (start, end)
    np.savez_compressed(outfile, **param)


def submit(loomfile, model, hapcode, chunk, outdir, email, queue, mem, walltime, systype, dryrun):
    LOG.warn('Loom file: %s' % loomfile)
    LOG.warn('Models: %s, %s' % (model[0], model[1]))
    LOG.warn('HPC system type: %s' % systype)
    if dryrun:
        LOG.warn('Showing submission script only')

    with loompy.connect(loomfile) as ds:
        gsurv = np.where(ds.ra.selected)[0]
        num_gsurv = len(gsurv)
        num_cells = ds.shape[1]
        LOG.warn('The number of selected genes: %d' % num_gsurv)
        LOG.warn('The number of selected cells: %d' % num_cells)
        LOG.warn('%d jobs will be submitted' % int(np.ceil(num_gsurv/chunk)))

    if systype == 'pbs':
        tot_layer = ''
        mat_layer = hapcode[0]
        processed = 0
        for idx_start in xrange(0, num_gsurv, chunk):
            idx_end = min(idx_start+chunk, num_gsurv-1)
            start = gsurv[idx_start]
            end = gsurv[idx_end]
            #if idx_end == num_gsurv-1:
            #    end += 1
            LOG.info('Start: %d, End %d' % (start, end))
            infile = os.path.join(outdir, '_chunk.%05d-%05d.npz' % (start, end))
            data_dict = dict()
            genes = gsurv[idx_start:idx_end]
            processed += len(genes)
            LOG.debug('Genes: %s' % ' '.join(genes.astype(str)))
            LOG.debug('Total %d genes submitted in this job' % len(genes))
            data_dict['shape'] = (len(genes), num_cells)
            with loompy.connect(loomfile) as ds:
                data_dict['GeneID'] = ds.ra.GeneID[genes]
                cur_chunk = dict()
                cur_chunk[tot_layer] = ds.layers[tot_layer][genes, :]
                cur_chunk[mat_layer] = ds.layers[mat_layer][genes, :]
                data_dict['counts'] = cur_chunk
                data_dict['size'] = ds.ca.size
                data_dict['selected'] = np.ones(len(genes))  # select all
                np.savez_compressed(infile, **data_dict)
            outfile = os.path.join(outdir, 'scase.%05d-%05d.param.npz' % (start, end))
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
        LOG.debug('Total %d genes were submitted' % processed)
        LOG.warn('Job submission complete')
    elif systype == 'pbs-loom':
        for idx_start in xrange(0, num_gsurv, chunk):
            idx_end = min(idx_start+chunk, num_gsurv-1)
            start = gsurv[idx_start]
            end = gsurv[idx_end]
            LOG.info('Start: %d, End %d' % (start, end))
            infile = os.path.join(outdir, '_chunk.%05d-%05d.loom' % (start, end))
            genes = gsurv[idx_start:idx_end]
            LOG.debug('Genes: %s' % ' '.join(genes.astype()))
            with loompy.connect(loomfile) as ds:
                with loompy.new(infile) as dsout:
                    for (ix, selection, view) in ds.scan(items=genes, axis=0):
                        LOG.debug('Genes in this view: %s' % ' '.join(selection.astype()))
                        dsout.add_columns(view.layers, col_attrs=view.col_attrs, row_attrs=view.row_attrs)
            outfile = os.path.join(outdir, 'scase.%05d-%05d.param.npz' % (start, end))
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
        LOG.warn('Job submission complete')
    elif 'lsf':
        raise NotImplementedError('LSF submission is not yet supported')
    else:
        raise RuntimeError('No plan to support other job scheduling system until we see many requests')


def collate(indir, loomfile, filetype, filename, model):
    if model[0] == 'null' and model[1] == 'null':
        raise RuntimeError('At least either of ASE or TGX model should be specified.')

    if filetype == "counts":
        # Get cell IDs
        LOG.info('Searching subdirectories for cells at %s' % indir)
        dlist = glob.glob(os.path.join(indir, '*/'))
        if len(dlist) > 0:  # Assuming indir/cellID/filename
            LOG.warn('%d subdirectories were found' % len(dlist))
            clist = [os.path.basename(d.rstrip('/')) for d in dlist]
            clist.sort()
            flist = [os.path.join(indir, c, filename) for c in clist]
            for f in flist:
                if not os.path.exists(f):
                    raise FileNotFoundError('%s does not exist. Consider providing a full file name' % f)
        else:  # If a subdirectory for each cell does not exist
            LOG.warn('No subdirectories were found')
            LOG.warn('Looking at %s directly for count files...' % indir)
            flist = glob.glob(os.path.join(indir, filename))  # filename should include wildcard in this case
            if len(flist) > 0:
                LOG.warn('%d files were found under %s' % (len(flist), indir))
                flist.sort()
                clist = [os.path.basename(f).split('.')[0] for f in flist]  # Assuming basename is cell ID
            else:
                raise FileNotFoundError('No files to collate')
        if len(clist) != len(flist):
            raise RuntimeError('The numbers of files and cells do not match')

        f = flist[0]
        with open(f) as fh:
            curline = fh.readline()
            item = curline.rstrip().split('\t')
            hapcodes = item[1:-1]
        LOG.warn('Haplotypes: %s' % '\t'.join(hapcodes))
        geneID = np.loadtxt(f, dtype=str, skiprows=1, usecols=0)
        LOG.warn('Number of genes: %d [%s %s ...]' % (len(geneID), geneID[0], geneID[1]))
        ds = loompy.new(loomfile)
        LOG.warn('A new loom file created: %s' % loomfile)
        LOG.warn('Populating loom file with TGX')
        for cix, f in enumerate(flist):
            new_column = np.loadtxt(f, skiprows=1, usecols=(-1,))
            cellID = clist[cix]
            ds.add_columns(np.matrix(new_column).T, row_attrs={'GeneID': geneID.astype(str)},
                           col_attrs={'CellID': np.array([cellID], dtype=str), 'size': np.array([new_column.sum()])})
            LOG.info('TGX loaded from %s' % f)
        LOG.warn('Populating loom file with ASE')
        for hix, h in enumerate(hapcodes):
            LOG.info('Loading ASE for Haplotype %s' % h)
            ds.layers[h] = 'float64'
            for cix, f in enumerate(flist):
                ds.layers[h][:, cix] = np.loadtxt(f, skiprows=1, usecols=(hix+1,))
                LOG.info('ASE loaded from %s' % f)
        ds.close()
        LOG.warn('Done. You may want to add more row_attrs or col_attrs to %s' % loomfile)

    elif filetype == 'params':
        flist = glob.glob(os.path.join(indir, '*.param.npz'))  # All the param files are assumed to be in indir
        LOG.info(os.path.join(indir, '*.param.npz'))  # filename is not used in collate --params
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
            LOG.warn('Loading %s' % f)
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
                # Add handling of new ASE models here!!
                else:
                    raise NotImplementedError('scBASE does not know how to process %s model results' % model[0])

                # Process TGX results
                if model[1] == 'null':
                    pass
                elif model[1] == 'pg':
                    LOG.info('Writing the TGX results by PG model')
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


