from setuptools import setup
import os
setup(name='zwatershed',
      version='0.2',
      description='Fast watersheds',
      url='https://github.com/TuragaLab/zwatershed',
      author='Chandan Singh',
      author_email='csinva@virginia.edu',
      license='MIT',
      install_requires=['cython','numpy','h5py'],	
      packages=['zwatershed'],
      zip_safe=False)
os.system('./make.sh')
