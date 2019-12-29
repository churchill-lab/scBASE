# -*- coding: utf-8 -*-

"""Top-level package for scBASE."""

import os
#from .scbase import *
from .scbase import load_model, run_em, run_mcmc, run_mcmc, __em4ase, __mcmc4ase, __em4tgx, __mcmc4tgx, submit, collate, adjust

__author__ = """Kwangbom "KB" Choi, Ph.D."""
__email__ = 'kb.choi@jax.org'
__version__ = '0.1.1'
__logo__ = """
┌─┐┌─┐╔╗ ╔═╗╔═╗╔═╗
└─┐│  ╠╩╗╠═╣╚═╗║╣ 
└─┘└─┘╚═╝╩ ╩╚═╝╚═╝
            v""" + __version__
