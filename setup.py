from setuptools import setup

setup(name='zwatershed',
      version='0.1',
      description='Fast watersheds',
      url='https://github.com/TuragaLab/zwatershed',
      author='Chandan Singh',
      author_email='csinva@virginia.edu',
      license='MIT',
      install_requires=['cython','numpy','h5py'],	
      packages=['zwatershed'],
      zip_safe=False)
