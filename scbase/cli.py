# -*- coding: utf-8 -*-

"""Console script for scbase."""
from __future__ import print_function
import sys
import click

from . import scbase
from . import utils
from . import __logo__, __version__
from . import get_data


@click.group()
@click.version_option(version=__version__, message=__logo__)
def main(args=None):
    pass


@main.command()
@click.argument('alnfile', metavar='<alnfile>', type=click.Path(exists=True, resolve_path=True, dir_okay=False))
@click.option('--start', metavar='<cix_start>', type=int, default=0, help='Starting cell (column index)')
@click.option('--end', metavar='<cix_end>', type=int, default=None, help='Ending cell (column index)')
@click.option('-v', '--verbose', count=True, help='the more times listed, the more output')
def disambiguate(alnfile, start, end, verbose):
    """
    Disambiguates multi-mapping reads
    :param alnfile:
    :param start:
    :param end:
    :param verbose:
    """
    utils.configure_logging(verbose)
    scbase.disambiguate(alnfile, start, end)


@main.command()
@click.argument('loomfile', metavar='<loomfile>', type=click.Path(exists=True, resolve_path=True, dir_okay=False))
@click.option('--min-read-count', metavar='<min_read_count>', type=int, default=0, help='Min read count required')
@click.option('--min-cell-count', metavar='<min_cell_count>', type=int, default=0, help='Min cell count required')
@click.option('--layer', metavar='<layer_key>', type=str, default='', help='The layer to consider selection')
@click.option('-v', '--verbose', count=True, help='the more times listed, the more output')
def select(loomfile, min_read_count, min_cell_count, layer, verbose):
    """
    Select genes if they are expressed at least <min_read_count> in at least <min_cell_count> cells
    :param loomfile:
    :param min_read_count:
    :param min_cell_count:
    :param layer:
    :param verbose:
    :return:
    """
    utils.configure_logging(verbose)
    scbase.select(loomfile, min_read_count, min_cell_count, layer)


@main.command()
@click.argument('loomfile', metavar='<loomfile>', type=click.Path(exists=True, resolve_path=True, dir_okay=False))
@click.option('--hapcode', metavar='<mat_hapcode> <pat_hapcode>', type=(str, str), default=('M', 'P'))
@click.option('-m', '--model', metavar='<ase_model> <tgx_model>', type=(str, str), default=('zoibb', 'pg'))
@click.option('-s', '--start', metavar='<gix_start>', type=int, default=0, help='Starting gene (row index)')
@click.option('-e', '--end', metavar='<gix_end>', type=int, default=None, help='Ending gene (row index)')
@click.option('-o', '--outfile', metavar='<outfile>', type=click.Path(resolve_path=True, dir_okay=False), default=None)
@click.option('-v', '--verbose', count=True, help='\'-v\' is Level 1 and \'-vv\' is Level 2')
def run_mcmc_from_loom(loomfile, model, hapcode, start, end, outfile, verbose):
    """
    Run MCMC (using loom file)
    :param loomfile:
    :param hapcode:
    :param model:
    :param start:
    :param end:
    :param outfile:
    :param verbose:
    """
    utils.configure_logging(verbose)
    scbase.run_mcmc_from_loom(loomfile, model, hapcode, start, end, outfile)


@main.command()
@click.argument('datafile', metavar='<npzfile>')
@click.option('--hapcode', metavar='<mat_hapcode> <pat_hapcode>', type=(str, str), default=('M', 'P'))
@click.option('-m', '--model', metavar='<ase_model> <tgx_model>', type=(str, str), default=('zoibb', 'pg'))
@click.option('-s', '--start', metavar='<gix_start>', type=int, default=0, help='Starting gene (row index)')
@click.option('-e', '--end', metavar='<gix_end>', type=int, default=None, help='Ending gene (row index)')
@click.option('-o', '--outfile', metavar='<outfile>', type=click.Path(resolve_path=True, dir_okay=False), default=None)
@click.option('-v', '--verbose', count=True, help='\'-v\' is Level 1 and \'-vv\' is Level 2')
def run_mcmc(datafile, model, hapcode, start, end, outfile, verbose):
    """
    Run MCMC (using npz file)
    :param datafile:
    :param hapcode:
    :param model:
    :param start:
    :param end:
    :param outfile:
    :param verbose:
    """
    utils.configure_logging(verbose)
    scbase.run_mcmc(datafile, model, hapcode, start, end, outfile)


@main.command()
@click.argument('loomfile', metavar='<loomfile>')
@click.option('--hapcode', metavar='<mat_hapcode> <pat_hapcode>', type=(str, str), default=('M', 'P'))
@click.option('-m', '--model', metavar='<ase_model> <tgx_model>', type=(str, str), default=('zoibb', 'pg'))
@click.option('-r', '--common-scale', metavar='<common_scale>', type=float, default=10000)
@click.option('-p', '--percentile', metavar='<percentile>', type=int, default=50)
@click.option('-s', '--start', metavar='<gix_start>', type=int, default=0, help='Starting gene (row index)')
@click.option('-e', '--end', metavar='<gix_end>', type=int, default=None, help='Ending gene (row index)')
@click.option('-t', '--tol', metavar='<tolerance>', type=float, default=0.000001)
@click.option('-i', '--max-iters', metavar='<max_iters>', type=int, default=100)
@click.option('-o', '--outfile', metavar='<outfile>', type=click.Path(resolve_path=True, dir_okay=False), default=None)
@click.option('-v', '--verbose', count=True, help='\'-v\' is Level 1 and \'-vv\' is Level 2')
def run_em(loomfile, model, common_scale, percentile, hapcode, start, end, tol, max_iters, outfile, verbose):
    """
    Run EM
    :param loomfile:
    :param model:
    :param common_scale:
    :param percentile:
    :param hapcode:
    :param start:
    :param end:
    :param tol:
    :param max_iters:
    :param outfile:
    :param verbose:
    :return:
    """
    utils.configure_logging(verbose)
    scbase.run_em(model, common_scale, percentile, hapcode, start, end, tol, max_iters, outfile)


@main.command()
@click.argument('loomfile', metavar='<loomfile>', type=click.Path(exists=True, dir_okay=False))
@click.option('--hapcode', metavar='<mat_hapcode> <pat_hapcode>', type=(str, str), default=('M', 'P'))
@click.option('-m', '--model', metavar='<ase_model> <tgx_model>', type=(str, str), default=('zoibb', 'pg'))
@click.option('-c', '--chunk', metavar='<chunk_size>', type=int, default=25, help='Chunk size')
@click.option('-o', '--outdir', metavar='<outdir>', type=click.Path(exists=True, resolve_path=True, file_okay=False), default='.')
@click.option('--systype', metavar='<systype>', default='pbs', help='Type of HPC cluster system (default: pbs)')
@click.option('--email', metavar='<email>', type=str, default=None, help='Notification E-mail')
@click.option('--queue', metavar='<queue>', type=str, default=None, help='Queue name')
@click.option('--mem', metavar='<mem>', type=int, default=0, help='Memory in GB (default: 16GB)')
@click.option('--walltime', metavar='<walltime>', type=int, default=0, help='Walltime in hours (default: 24h')
@click.option('--dryrun', is_flag=True)
@click.option('-v', '--verbose', count=True, help='\'-v\' is Level 1 and \'-vv\' is Level 2')
def submit(loomfile, model, hapcode, chunk, outdir, email, queue, mem, walltime, systype, dryrun, verbose):
    """
    Submit jobs to HPC clusters
    :param loomfile:
    :param model:
    :param hapcode:
    :param chunk:
    :param outdir:
    :param email:
    :param queue:
    :param mem:
    :param walltime:
    :param systype:
    :param dryrun:
    :param verbose:
    :return:
    """
    utils.configure_logging(verbose)
    scbase.submit(loomfile, model, hapcode, chunk, outdir, email, queue, mem, walltime, systype, dryrun)


@main.command()
@click.argument('indir', metavar='<indir>', type=click.Path(file_okay=False))
@click.argument('loomfile', metavar='<loomfile>', type=click.Path(dir_okay=False))
@click.option('--counts', 'filetype', flag_value='counts')
@click.option('--params', 'filetype', flag_value='params')
@click.option('--name', 'filename', default='*genes*counts')
@click.option('-m', '--model', metavar='<ase_model> <tgx_model>', type=(str, str), default=('zoibb', 'pg'))
@click.option('-v', '--verbose', count=True, help='\'-v\' is Level 1 and \'-vv\' is Level 2')
def collate(loomfile, indir, filetype, filename, model, verbose):
    """
    Utility function that collates input/output files
    :param loomfile:
    :param indir:
    :param filetype:
    :param filename:
    :param model:
    :param verbose:
    :return:
    """
    utils.configure_logging(verbose)
    scbase.collate(indir, loomfile, filetype, filename, model)


if __name__ == "__main__":
    sys.exit(main())
