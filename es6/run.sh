#!/bin/bash

cd build
make -j4
cd ..
./build/prng/main 
