"""
scQTLib: single-cell quantitative trait analysis library
See: https://github.com/StatBiomed/scTSS
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Set __version__ for the project.
exec(open("./scTSS/version.py").read())

# Get the long description from the relevant file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()
    
reqs = ['numpy>=1.9.0', 'scipy>=1.4.0', 'matplotlib','anndata>=0.6',
'scanpy>=1.5','pysam>=0.15.2','brie>=2.2.0','pandas>=0.23.0','scikit-learn>=0.23']

setup(
    name='scTSS',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=__version__,

    description='scTSS: Detection and couting alternative TSS in single cells',
    long_description=long_description,

    # The project's main homepage.
    url='https://github.com/StatBiomed/scTSS',

    # Author details
    author=['scTSS Team'],
    author_email='',

    # Choose your license
    license='Apache-2.0',

    # What does your project relate to?
    keywords=['Transcript Start Site', 'single-cell RNA-seq'],

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(),

    entry_points={
          'console_scripts': [
            'scTSS = scTSS.bin.scTSS_main:main',
            'scTSS-count=scTSS.bin.count:main',
            'scTSS-quant=scTSS.bin.quant:main',
            'scTSS-plot=scTSS.bin.plot:main'
            ],
          }, 

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html

    python_requires='>=3.5',

    install_requires=reqs,

    extras_require={
        'docs': [
            #'sphinx == 1.8.3',
            'sphinx_bootstrap_theme']},

    py_modules = ['scTSS']

    # buid the distribution: python setup.py sdist
    # upload to pypi: twine upload dist/...
)