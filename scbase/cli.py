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
@click.argument('alntools_file', metavar='alnfile', type=click.Path(exists=True, resolve_path=True, dir_okay=False))
@click.option('--start', metavar='c_start', default=0, help='Starting cell (column index)')
@click.option('--end', metavar='c_end', default=None, help='Ending cell (column index)')
@click.option('-v', '--verbose', count=True, help='the more times listed, the more output')
def disambiguate(alnfile, c_start, c_end, verbose):
    """Console script for scbase
    :param alnfile:
    :param c_start:
    :param c_end:
    :param verbose:
    """
    utils.configure_logging(verbose)
    scbase.disambiguate(alnfile, c_start, c_end)


@main.command()
@click.argument('loomfile', metavar='loomfile', type=click.Path(exists=True, resolve_path=True, dir_okay=False))
@click.option('-m', '--model', metavar='model', default='zoibb')
@click.option('--start', metavar='g_start', default=0, help='Starting gene (row index)')
@click.option('--end', metavar='g_end', default=None, help='Ending gene (row index)')
@click.option('-v', '--verbose', count=True, help='\'-v\' is Level 1 and \'-vv\' is Level 2')
def run_mcmc(loomfile, model, g_start, g_end, verbose):
    """MCMC script for scBASE
    :param loomfile:
    :param model:
    :param g_start:
    :param g_end:
    :param verbose:
    """
    utils.configure_logging(verbose)
    scbase.run_mcmc(loomfile, model, g_start, g_end)


@main.command()
@click.argument('indir', metavar='indir', type=click.Path(exists=True, resolve_path=True))
@click.argument('loomfile', metavar='loomfile', type=click.Path(exists=True, resolve_path=True, dir_okay=False))
@click.option('--counts', 'filetype', flag_value='counts')
@click.option('--params', 'filetype', flag_value='params')
@click.option('-n', '--name', metavar='filename', default=None)
@click.option('-v', '--verbose', count=True, help='\'-v\' is Level 1 and \'-vv\' is Level 2')
def collate(loomfile, indir, filetype, filename, verbose):
    """
    Utility function that collates input/output files
    :param loomfile:
    :param indir:
    :param filetype:
    :param filename:
    :param verbose:
    :return:
    """
    utils.configure_logging(verbose)
    scbase.collate(indir, loomfile, filetype, filename)


if __name__ == "__main__":
    sys.exit(main())
