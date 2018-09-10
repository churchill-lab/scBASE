# -*- coding: utf-8 -*-

"""Main module."""

from . import utils
from . import get_data

LOG = utils.get_logger()

try:
    xrange
except NameError:
    xrange = range


def disambiguate(alntools_file):
    LOG.warn('Quantifying allele-specific expression in each cell')
    LOG.info('Level-1 verbose')
    LOG.debug('Level-2 verbose')
    print(get_data('README'))

