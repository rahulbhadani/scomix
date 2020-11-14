import setuptools
from distutils.dir_util import copy_tree
from pathlib import Path
PACKAGE_NAME='scomix'
import shutil, os
shutil.copy('README.md', PACKAGE_NAME + '/README.md')

def readme():
    with open("README.md", "r") as fh:
        long_description = fh.read()
        return long_description

v = Path(PACKAGE_NAME + "/version").open(encoding = "utf-8").read().splitlines()


setuptools.setup(
    name='scomix',
    version=v[0].strip(),
    author="Rahul Bhadani",
    author_email="rahulbhadani@email.arizona.edu",
    description="Python package for analysis of data from single-cell data-science.",
    long_description=readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/rahulbhadani/scomix",
    packages=setuptools.find_packages(),
    install_requires=[
        l.strip() for l in Path("requirements.txt").open(encoding = "utf-8").read().splitlines()
        ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Framework :: IPython",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "License :: OSI Approved :: MIT License",
        ],
    keywords='single-cell, scRNAseq, transcriptomics, atacseq, rnaseq, bioinformatics, life-sciences, computation-biology, data-science, system-biology',
    include_package_data=True,
    package_data={'scomix': ['README.md','version']},
    zip_safe=False
        )

os.remove('scomix/README.md')

