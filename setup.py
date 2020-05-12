from setuptools import setup
import os
import sys

project_dir = os.path.abspath(os.path.dirname(__file__))

if sys.version_info[0] < 3:
    with open(os.path.join(project_dir, 'README.rst')) as f:
        long_description = f.read()
else:
    with open(os.path.join(project_dir, 'README.rst'), encoding='utf-8') as f:
        long_description = f.read()

setup(
    name='gxfutils',
	version='1.0.0',    
    description=('A collection of utilities for manipulating GTF/GFF files.'),
    long_description=long_description,
    author='Bimber Lab',
    author_email='bimber@ohsu.edu',
    url='https://github.com/bimberlabinternal/GxfUtils',
    license='MPL-2.0',
	install_requires=project_dir.joinpath('requirements.txt').read_text().split('\n'),
    include_package_data=True,
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6'
	],    
	entry_points = {
        'console_scripts': [
            'gxfutils = gxfutils.__main__:main'
        ]
    }
)