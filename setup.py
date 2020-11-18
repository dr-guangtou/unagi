# -*- coding: utf-8 -*-

'''This sets up the package.
Stolen from http://python-packaging.readthedocs.io/en/latest/everything.html and
modified by me.
'''
from setuptools import setup, find_packages

def readme():
    """Load the README file."""
    with open('README.md') as f:
        return f.read()

# let's be lazy and put requirements in one place
# what could possibly go wrong?
with open('requirements.txt') as infd:
    INSTALL_REQUIRES = [x.strip('\n') for x in infd.readlines()]


###############
## RUN SETUP ##
###############

setup(
    name='unagi',
    description=('Delicious Data from Hyper Suprime-Cam Survey!'),
    long_description=readme(),
    long_description_content_type="text/markdown",
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: MIT License',
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ],
    keywords='astronomy',
    url='https://github.com/dr-guangtou/unagi',
    author='Song Huang',
    author_email='shuang89@ucsc.edu',
    license='MIT',
    packages=find_packages(),
    package_data={
        'unagi': ["data/camera/*", "data/filters/*",
                  "data/pdr1_wide/*", "data/pdr2_dud/*", "data/pdr2_wide/*",
                  "data/solar/*", "data/sql_template/*"]
    },
    install_requires=INSTALL_REQUIRES,
    include_package_data=True,
    zip_safe=False,
    python_requires='>=3.6',
    scripts=['bin/hsc_bulk_cutout'],
    use_scm_version = {
        "local_scheme": "no-local-version"
    },
    setup_requires=['setuptools_scm'],
)
