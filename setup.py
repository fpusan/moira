# -*- coding: utf-8 -*-
from distutils.core import setup, Extension



bernoulli = Extension('bernoulli', sources = ['moira/bernoullimodule.c'])
nw_align = Extension('nw_align', sources = ['moira/nw_align.c'])

setup(
    name = 'moira',
    packages = ['moira'],
    version = 'v1.0.0',
    description = 'A random test lib',
    author = 'Fernando Puente-Sánchez',
    author_email = 'fpusan@gmail.com',
    url = 'https://github.com/fpusan/moira', # use the URL to the github repo
    download_url = 'https://github.com/fpusan/moira/tarball/1.0.0', # I'll explain this in a second
    license='BSD',
    keywords = ['testing', 'logging', 'example'], # arbitrary keywords
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',

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
)
