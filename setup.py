#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

import os
import scbase
from setuptools import setup, find_packages
try:
    from setuptools.command.install import install
    from setuptools.command.build_py import build_py
except ImportError:
    from distutils.command.install import install
    from distutils.command.build_py import build_py

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['numpy',
                'scipy',
                'Click>=6.0',
                'pystan>=2.17',
                'h5py>=2.8',
                'pandas>=0.23.4',
                'loompy>=2.0.14',
                'future',
                'six']  # 'alntools>=0.1.0',
setup_requirements = [ ]
test_requirements = [ ]


class StanBuild(build_py):
    def run(self):
        build_py.run(self)
        build_path = os.path.join(self.build_lib, 'scbase/stan/')
        self.mkpath(build_path)
        target_files = []

        def compile():
            import glob
            import pystan
            try:
                import cPickle as pickle
            except:
                import pickle
            stan_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'scbase/stan')
            flist = glob.glob(os.path.join(stan_path, '*.stan'))
            compile_flags = [
                '-O3',
                '-ftemplate-depth-256',
                '-Wno-unused-function',
                '-Wno-uninitialized',
            ]
            for f in flist:
                fbase = os.path.splitext(f)[0]
                fout = fbase + ".pkl"
                target_files.append(fout)
                if os.path.exists(fout):
                    print("Compiled stan file exists: %s" % fout)
                else:
                    print("Compiling %s..." % f)
                    stan_model = pystan.StanModel(file=f, extra_compile_args=compile_flags)
                    with open(fout, 'wb') as fhout:
                        pickle.dump(stan_model, fhout)

        self.execute(compile, [], 'Compiling STAN code:')

        # copy resulting tool to library build folder
        print(build_path)
        if not self.dry_run:
            for tfile in target_files:
                self.copy_file(tfile, build_path)


setup(
    author="Kwangbom \"KB\" Choi, Ph.D.",
    author_email='kb.choi@jax.org',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description="**scBASE** is a python implementation of \"soft\" zero-and-one inflated model for estimating cellular allelic proportions from scRNA-Seq data",
    entry_points={
        'console_scripts': [
            'scbase=scbase.cli:main',
        ],
    },
    install_requires=requirements,
    cmdclass={'build_py': StanBuild},
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='scbase',
    name='scbase',
    packages=find_packages(include=['scbase']),
    setup_requires=setup_requirements,
    scripts=[
        'scripts/run_mcmc_on_cluster.sh',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/churchill-lab/scBASE',
    version=scbase.__version__,
    zip_safe=False,
)
