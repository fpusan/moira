# -*- coding: utf-8 -*-
from distutils.core import setup, Extension, Command
import sys, subprocess

class TestMoira(Command):
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        errno = subprocess.call([sys.executable, '-m', 'test.test_moira'], cwd = 'moira')
        raise SystemExit(errno)


bernoulli = Extension('bernoulli', sources = ['moira/bernoullimodule.c'])
nw_align = Extension('nw_align', sources = ['moira/nw_align.c'])

setup(
    name = 'moira',
    packages = ['moira'],
    version = 'v1.3.0',
    description = 'Quality-filter raw sequence reads using the Poisson binomial filtering algorithm',
    author = 'Fernando Puente-Sánchez',
    author_email = 'fpusan@gmail.com',
    url = 'https://github.com/fpusan/moira', # use the URL to the github repo
    download_url = 'https://github.com/fpusan/moira/tarball/v1.3.0',
    license = 'BSD-3',
    keywords = ['high-throughput sequencing', 'microbial ecology', '16S analysis', 'marker-gene',
                'bioinformatics', 'Illumina', '454', 'IonTorrent', 'quality-filtering'],
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 5 - Production/Stable',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        # Pick your license as you wish (should match "license" above)
         'License :: OSI Approved :: BSD License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
    ],
    ext_modules = [bernoulli, nw_align],
    scripts = ['moira/moira.py'],
    cmdclass = {'test': TestMoira}
)
