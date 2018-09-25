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
    """Console script for scbase
    :param alnfile:
    :param start:
    :param end:
    :param verbose:
    """
    utils.configure_logging(verbose)
    scbase.disambiguate(alnfile, start, end)


@main.command()
@click.argument('loomfile', metavar='<loomfile>', type=click.Path(exists=True, resolve_path=True, dir_okay=False))
@click.option('--hapcode', metavar='<mat_hapcode> <pat_hapcode>', type=(str, str), default=('M', 'P'))
@click.option('-m', '--model', metavar='<ase_model> <tgx_model>', type=(str, str), default=('zoibb', 'pg'))
@click.option('-s', '--start', metavar='<gix_start>', type=int, default=0, help='Starting gene (row index)')
@click.option('-e', '--end', metavar='<gix_end>', type=int, default=None, help='Ending gene (row index)')
@click.option('-o', '--outfile', metavar='<outfile>', type=click.Path(resolve_path=True, dir_okay=False), default=None)
@click.option('-v', '--verbose', count=True, help='\'-v\' is Level 1 and \'-vv\' is Level 2')
def run_mcmc(loomfile, model, hapcode, start, end, outfile, verbose):
    """MCMC script for scBASE
    :param loomfile:
    :param hapcode:
    :param model:
    :param start:
    :param end:
    :param outfile:
    :param verbose:
    """
    utils.configure_logging(verbose)
    scbase.run_mcmc(loomfile, model, hapcode, start, end, outfile)


@main.command()
@click.argument('loomfile', metavar='<loomfile>', type=click.Path(dir_okay=False))
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
