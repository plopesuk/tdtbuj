#!/bin/bash

cd ../
rm -rf CMakeFiles
rm -f CMakeCache.txt 
rm -f cmake_install.cmake
rm -f Makefile
rm -rf bin/*
cmake .
make
cd src
