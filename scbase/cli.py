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


# @main.command()
# @click.argument('alnfile', metavar='<alnfile>', type=click.Path(exists=True, resolve_path=True, dir_okay=False))
# @click.option('--start', metavar='<cix_start>', type=int, default=0, help='Starting cell (column index)')
# @click.option('--end', metavar='<cix_end>', type=int, default=None, help='Ending cell (column index)')
# @click.option('-v', '--verbose', count=True, help='the more times listed, the more output')
# def disambiguate(alnfile, start, end, verbose):
#     """
#     Disambiguates multi-mapping reads
#     """
#     utils.configure_logging(verbose)
#     scbase.disambiguate(alnfile, start, end)


@main.command()
@click.argument('loomfile', metavar='<loomfile>', type=click.Path(exists=True, resolve_path=True, dir_okay=False))
@click.option('--min-read-count', metavar='<min_read_count>', type=int, default=0, help='Min read count required')
@click.option('--min-cell-count', metavar='<min_cell_count>', type=int, default=0, help='Min cell count required')
@click.option('--layer', metavar='<layer_key>', type=str, default='', help='The layer to consider selection')
@click.option('-v', '--verbose', count=True, help='\'-v\' is Level 1 and \'-vv\' is Level 2')
def select(loomfile, min_read_count, min_cell_count, layer, verbose):
    """
    Selects genes if they are expressed at least <min_read_count> in at least <min_cell_count> cells
    """
    utils.configure_logging(verbose)
    scbase.select(loomfile, min_read_count, min_cell_count, layer)


@main.command()
@click.argument('loomfile', metavar='<loomfile>', type=click.Path(exists=True, resolve_path=True, dir_okay=False))
@click.option('--hapcode', metavar='<mat_hapcode> <pat_hapcode>', type=(str, str), default=('M', 'P'),
help='Haplotype code for maternal and paternal alleles (default: M P)')
@click.option('-m', '--model', metavar='<ase_model> <tgx_model>', type=(str, str), default=('zoibb', 'pg'),
help='This is a developer option used only when other models in addition to ones provided by default are available.')
@click.option('-s', '--start', metavar='<gix_start>', type=int, default=0, help='Starting gene (row index)')
@click.option('-e', '--end', metavar='<gix_end>', type=int, default=None, help='Ending gene (row index)')
@click.option('-o', '--outfile', metavar='<outfile>', type=click.Path(resolve_path=True, dir_okay=False), default=None,
help='Name of output file (default: scbase.<start>-<end>.param.npz)')
@click.option('-v', '--verbose', count=True, help='\'-v\' is Level 1 and \'-vv\' is Level 2')
def run_mcmc(loomfile, model, hapcode, start, end, outfile, verbose):
    """
    Runs MCMC (using loom file)
    """
    utils.configure_logging(verbose)
    scbase.run_mcmc(loomfile, model, hapcode, start, end, outfile)


@main.command()
@click.argument('datafile', metavar='<npzfile>')
@click.option('--hapcode', metavar='<mat_hapcode> <pat_hapcode>', type=(str, str), default=('M', 'P'),
help='Haplotype code for maternal and paternal alleles (default: M P)')
@click.option('-m', '--model', metavar='<ase_model> <tgx_model>', type=(str, str), default=('zoibb', 'pg'),
help='This is a developer option used only when other models in addition to ones provided by default are available.')
@click.option('-s', '--start', metavar='<gix_start>', type=int, default=0, help='Starting gene (row index)')
@click.option('-e', '--end', metavar='<gix_end>', type=int, default=None, help='Ending gene (row index)')
@click.option('-o', '--outfile', metavar='<outfile>', type=click.Path(resolve_path=True, dir_okay=False), default=None,
help='Name of output file (default: scbase.<start>-<end>.param.npz)')
@click.option('-v', '--verbose', count=True, help='\'-v\' is Level 1 and \'-vv\' is Level 2')
def run_mcmc_from_npz(datafile, model, hapcode, start, end, outfile, verbose):
    """
    Runs MCMC (using npz file)
    """
    utils.configure_logging(verbose)
    scbase.run_mcmc_from_npz(datafile, model, hapcode, start, end, outfile)


@main.command()
@click.argument('loomfile', metavar='<loomfile>')
@click.option('--hapcode', metavar='<mat_hapcode> <pat_hapcode>', type=(str, str), default=('M', 'P'),
help='Haplotype code for maternal and paternal alleles (default: M P)')
@click.option('-m', '--model', metavar='<ase_model> <tgx_model>', type=(str, str), default=('zoibb', 'pg'),
help='This is a developer option used only when other models in addition to ones provided by default are available.')
@click.option('-r', '--common-scale', metavar='<common_scale>', type=float, default=10000,
help='Read counts per cell after scaliing (default: 10000)')
@click.option('-p', '--percentile', metavar='<percentile>', type=int, default=50)
@click.option('-s', '--start', metavar='<gix_start>', type=int, default=0, help='Starting gene (row index)')
@click.option('-e', '--end', metavar='<gix_end>', type=int, default=None, help='Ending gene (row index)')
@click.option('-t', '--tol', metavar='<tolerance>', type=float, default=0.000001, help='Tolerance for termination (default: 0.000001)')
@click.option('-i', '--max-iters', metavar='<max_iters>', type=int, default=100, help='Max iterations for termination (default: 100)')
@click.option('-o', '--outfile', metavar='<outfile>', type=click.Path(dir_okay=False), default=None, help='Name of outpput file')
@click.option('-v', '--verbose', count=True, help='\'-v\' is Level 1 and \'-vv\' is Level 2')
def run_em(loomfile, model, common_scale, percentile, hapcode, start, end, tol, max_iters, outfile, verbose):
    """
    Runs EM
    """
    utils.configure_logging(verbose)
    scbase.run_em(loomfile, model, common_scale, percentile, hapcode, start, end, tol, max_iters, outfile)


@main.command()
@click.argument('loomfile', metavar='<loomfile>', type=click.Path(exists=True, dir_okay=False))
@click.option('--hapcode', metavar='<mat_hapcode> <pat_hapcode>', type=(str, str), default=('M', 'P'),
help='Haplotype code for maternal and paternal alleles (default: M P)')
@click.option('-m', '--model', metavar='<ase_model> <tgx_model>', type=(str, str), default=('zoibb', 'pg'),
help='This is a developer option used only when other models in addition to ones provided by default are available.')
@click.option('-c', '--chunk', metavar='<chunk_size>', type=int, default=25, help='Number of genes in each chunk')
@click.option('-o', '--outdir', metavar='<outdir>', type=click.Path(exists=True, resolve_path=True, file_okay=False), default='.',
help='Folder name to store parameter files')
@click.option('--systype', metavar='<systype>', default='pbs', help='Type of HPC cluster system (default: pbs)')
@click.option('--email', metavar='<email>', type=str, default=None, help='Notification E-mail')
@click.option('--queue', metavar='<queue>', type=str, default=None, help='Queue name')
@click.option('--mem', metavar='<mem>', type=int, default=0, help='Memory in GB (default: 16GB)')
@click.option('--walltime', metavar='<walltime>', type=int, default=0, help='Walltime in hours (default: 24h)')
@click.option('--dryrun', is_flag=True, help='Use this when you want to rehearse your submit commands')
@click.option('-v', '--verbose', count=True, help='\'-v\' is Level 1 and \'-vv\' is Level 2')
def submit(loomfile, model, hapcode, chunk, outdir, email, queue, mem, walltime, systype, dryrun, verbose):
    """
    Submits scBASE fitting jobs to HPC clusters
    """
    utils.configure_logging(verbose)
    scbase.submit(loomfile, model, hapcode, chunk, outdir, email, queue, mem, walltime, systype, dryrun)


@main.command()
@click.argument('indir', metavar='<indir>', type=click.Path(file_okay=False))
@click.argument('loomfile', metavar='<loomfile>', type=click.Path(dir_okay=False))
@click.option('--counts', 'filetype', flag_value='counts', help='If you are collating count files')
@click.option('--params', 'filetype', flag_value='params', help='If you are collating param files')
@click.option('--name', 'filename', help='File name (Enclose the globbing characters with quotes)')
@click.option('-t', '--tidfile', metavar='<tidfile>', type=click.Path(resolve_path=True, dir_okay=False), default=None,
help='Name of target ID file')
@click.option('-m', '--model', metavar='<ase_model> <tgx_model>', type=(str, str), default=('zoibb', 'pg'),
help='This is a developer option used only when other models in addition to ones provided by default are available.')
@click.option('-v', '--verbose', count=True, help='\'-v\' is Level 1 and \'-vv\' is Level 2')
def collate(indir, loomfile, tidfile, filetype, filename, model, verbose):
    """
    Collates count or parameter files and store in a loom file (http://loompy.org).
    """
    utils.configure_logging(verbose)
    scbase.collate(indir, loomfile, tidfile, filetype, filename, model)


if __name__ == "__main__":
    sys.exit(main())
