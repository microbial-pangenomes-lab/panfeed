from setuptools import setup, find_packages
from codecs import open
from os import path
import os
import re
import io


def read(*names, **kwargs):
    with io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")
    ) as fp:
        return fp.read()


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


here = path.abspath(path.dirname(__file__))


with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='panfeed',
    version=find_version("panfeed/__init__.py"),
    description='Compute gene-cluster specific k-mers over a pangenome',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/microbial-pangenomes-lab/panfeed',
    author='Hannes Neubauer',
    author_email='Neubauer.Hannes@mh-hannover.de',
    license='Apache',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    packages=['panfeed'],
    entry_points={
        "console_scripts": [
            'panfeed = panfeed.__main__:main',
            'panfeed-get-clusters = panfeed.get_clusters:main',
            'panfeed-get-kmers = panfeed.get_kmers:main',
            'panfeed-plot = panfeed.plot:main',
            ]
    },
    install_requires=['numpy',
                      'pandas',
                      'pyfaidx',
                      'matplotlib',
                      'seaborn']
    #test_suite="tests",
)
