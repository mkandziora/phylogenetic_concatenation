#!/usr/bin/env python

from setuptools import setup

setup(name='concat',
      version='0.1',
      description='Concatenation of phylogenetic data',
      author='Martha Kandziora',
      author_email='martha.kandziora@yahoo.com',
      packages=['Concat'],
      install_requires=[
          'dendropy'
      ]
     )
