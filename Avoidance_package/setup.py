#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 10:25:44 2019

@author: bikash
"""

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()



setuptools.setup(
        name="avoidance2",
        version="0.0.0.9.1",
        author="Bikash<bikash.bhandari@postgrad.otago.ac.nz>, Lim<chunshen.lim@otago.ac.nz>",
        author_email="mk.bikash@gmail.com",
        description="Core algorithms for RNA sampling and optimization",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://github.com/Gardner-BinfLab/Avoidance-v2",
        packages=setuptools.find_packages(),
        include_package_data=True,
        classifiers=[
                "Development Status :: 4 - Beta",
                "Programming Language :: Python :: 3.6",
                'Intended Audience :: Science/Research',
                "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
                "Operating System :: POSIX :: Linux",
                'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
        keywords="RNA, MarkovModel, Optimization",
        install_requires=[
                'pandas==0.23.4',
                'scikit-learn==0.20.2'
                ],
        python_requires='>=3.6',
        zip_safe=False
)
