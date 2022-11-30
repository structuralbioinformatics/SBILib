from setuptools import setup, find_packages  # Always prefer setuptools over distutils
import os
import SBI

version = SBI.__version__


def read_file(path):
    with open(os.path.join(os.path.dirname(__file__), path)) as fp:
        return fp.read()


setup(
    name='StrBioInfo',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # http://packaging.python.org/en/latest/tutorial.html#version
    version=version,

    description='The StructuralBioInformatics Library',
    long_description=read_file('README.rst'),


    # Author details
    author='Jaume Bonet',
    author_email='jaume.bonet@gmail.com',

    # Choose your license
    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Operating System :: MacOS',
        'Operating System :: Unix',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],

    platforms='UNIX',
    keywords=['structural biology', 'bioinformatics'],

    install_requires=[x.strip() for x in open('REQUIREMENTS').readlines()],

    packages=find_packages(exclude=['docs', 'test']),
    include_package_data=True,

    package_data={
        'SBI': ['REQUIREMENTS'],
    },

    zip_safe=False,
)
