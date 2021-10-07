from setuptools import find_packages, setup

setup(name='fastatools',
      version='0.1.0',
      description='Fasta processing functions',
      url='#',
      author='Elliot Hill',
      install_requires=['Bio', 'pandas'],
      author_email='',
      packages=find_packages(include=['fastatools']),
      zip_safe=False)
