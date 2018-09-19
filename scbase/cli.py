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
@click.argument('alntools_file', metavar='<alnfile>', type=click.Path(exists=True, resolve_path=True, dir_okay=False))
@click.option('--start', metavar='<c_start>', default=0, help='Starting cell (column index)')
@click.option('--end', metavar='<c_end>', default=None, help='Ending cell (column index)')
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
@click.option('-m', '--model', metavar='<model>', type=(str, str), default=('zoibb', 'pg'))
@click.option('-s', '--start', metavar='<g_start>', default=0, help='Starting gene (row index)')
@click.option('-e', '--end', metavar='<g_end>', default=None, help='Ending gene (row index)')
@click.option('--maternal', metavar='<maternal>', default='M', help='Letter code for maternal haplotype')
@click.option('--paternal', metavar='<maternal>', default='M', help='Letter code for maternal haplotype')
@click.option('-v', '--verbose', count=True, help='\'-v\' is Level 1 and \'-vv\' is Level 2')
def run_mcmc(loomfile, model, start, end, maternal, verbose):
    """MCMC script for scBASE
    :param loomfile:
    :param model:
    :param maternal:
    :param start:
    :param end:
    :param verbose:
    """
    utils.configure_logging(verbose)
    scbase.run_mcmc(loomfile, model, maternal, start, end)


@main.command()
@click.argument('indir', metavar='<indir>', type=click.Path(file_okay=False))
@click.argument('loomfile', metavar='<loomfile>', type=click.Path(dir_okay=False))
@click.option('--counts', 'filetype', flag_value='counts')
@click.option('--params', 'filetype', flag_value='params')
@click.option('--name', 'filename', default=None)
@click.option('--group', 'groupname', default='group1')
@click.option('-m', '--model', metavar='<model>', type=(str, str), default=('zoibb', 'pg'))
@click.option('-v', '--verbose', count=True, help='\'-v\' is Level 1 and \'-vv\' is Level 2')
def collate(loomfile, indir, filetype, filename, groupname, model, verbose):
    """
    Utility function that collates input/output files
    :param loomfile:
    :param indir:
    :param filetype:
    :param filename:
    :param groupname:
    :param model:
    :param verbose:
    :return:
    """
    utils.configure_logging(verbose)
    scbase.collate(indir, loomfile, filetype, filename, groupname, model)


if __name__ == "__main__":
    sys.exit(main())
