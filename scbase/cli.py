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
@click.argument('alntools_file', metavar='alntools_file', type=click.Path(exists=True, resolve_path=True, dir_okay=False))
@click.option('-v', '--verbose', count=True, help='the more times listed, the more output')
def disambiguate(alntools_file, verbose):
    """Console script for scbase.
    :param alntools_file:
    :param verbose:
    """
    utils.configure_logging(verbose)
    scbase.disambiguate(alntools_file)


if __name__ == "__main__":
    sys.exit(main())
