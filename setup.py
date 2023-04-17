"""
CamoTSS: single-cell alternative TSS detection and analysis library
See: https://github.com/StatBiomed/CamoTSS
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path



here = path.abspath(path.dirname(__file__))

# Set __version__ for the project.
exec(open("./CamoTSS/version.py").read())

# Get the long description from the relevant file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()
    
reqs = ['numpy>=1.9.0', 'scipy>=1.4.0', 'matplotlib','anndata>=0.6','pyranges>=0.0.115',
'scanpy>=1.5','pysam>=0.15.2','brie>=2.2.0','pandas>=0.23.0','scikit-learn>=0.23','editdistance>=0.3.1']

setup(
    name='CamoTSS',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=__version__,

    description='CamoTSS: Detection alternative TSS in single cells',
    long_description=long_description,

    # The project's main homepage.
    url='https://github.com/StatBiomed/CamoTSS',

    # Author details
    author=['Ruiyan Hou'],
    author_email='ruiyan@connect.hku.hk',

    # Choose your license
    license='Apache-2.0',

    # What does your project relate to?
    keywords=['Transcript Start Site', 'single-cell RNA-seq'],

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().




    packages=find_packages(),

    package_data={'CamoTSS': ['model/*.sav']},

    entry_points={
          'console_scripts': [
            'CamoTSS = CamoTSS.bin.count:main',
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

    py_modules = ['CamoTSS']

    # buid the distribution: python setup.py sdist
    # upload to pypi: twine upload dist/...
)
