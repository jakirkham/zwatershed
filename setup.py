from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext
import os

class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)
        # Prevent numpy from thinking it is still in its setup process:
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())

setup(name='zwatershed',
      version='0.5',
      description='Fast watersheds',
      url='https://github.com/TuragaLab/zwatershed',
      author='Chandan Singh',
      author_email='csinva@virginia.edu',
	  cmdclass=dict(
            build_ext=build_ext
    	),
      license='MIT',
      install_requires=['cython','numpy','h5py'],	
	  setup_requires=['cython','numpy'],	
      packages=['zwatershed'],
	  ext_modules = [Extension("zwatershed/zwatershed",
                         ["zwatershed/zwatershed.pyx", "zwatershed/zwatershed_main.cpp"],
                         library_dirs=["zwatershed"],
                         language='c++',
                         # std= 'c++11',
                         extra_link_args=["-std=c++11"],
                         extra_compile_args=["-std=c++11", "-w"])],
      zip_safe=False)
