#!/usr/bin/env bash
cd zwatershed
rm zwatershed.so
python setup.py build_ext --inplace
printf "BUILD COMPLETE\n"
cd ..
