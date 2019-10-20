#!/bin/bash

#LGC test
c++ ./prng/die_lgc.cpp ./prng/prng.cpp -o ./prng/die_lgc
./prng/die_lgc |dieharder -a -g 200

#xorshiro128p test
c++ ./prng/die_xorshiro128p.cpp ./prng/prng.cpp ./prng/xorshiro.cpp -o ./prng/die_xorshiro128p
./prng/die_xorshiro128p |dieharder -a -g 200