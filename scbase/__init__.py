# -*- coding: utf-8 -*-

"""Top-level package for scBASE."""

import os


__author__ = """Kwangbom "KB" Choi, Ph.D."""
__email__ = 'kb.choi@jax.org'
__version__ = '0.1.0'
__logo__ = """
         _____ _____ _____ _____
 ___ ___| __  |  _  |   __|   __|
|_ -|  _| __ -|     |__   |   __|
|___|___|_____|__|__|_____|_____|
                           v""" + __version__

_ROOT = os.path.abspath(os.path.dirname(__file__))


def get_data(path):
    return os.path.join(_ROOT, 'stan', path)

