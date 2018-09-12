# -*- coding: utf-8 -*-

"""Main module."""

from . import utils
from . import get_data
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


def disambiguate(alntools_file):
    LOG.warn('Quantifying allele-specific expression in each cell')
    LOG.info('Level-1 verbose')
    LOG.debug('Level-2 verbose')
    print(get_data('README'))


def run_mcmc(alntools_file, model):
    LOG.warn('Quantifying allele-specific expression in each cell')
    LOG.info('Level-1 verbose')
    LOG.debug('Level-2 verbose')
    model_file = '%s.pkl' % model
    print(get_data(model_file))
    stan_model = pickle.load(open(get_data(model_file), 'rb'))
    print(stan_model.model_code)
