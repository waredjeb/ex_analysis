#!/bin/bash

#LGC test
./prng/die_lgc |dieharder -a -g 200

#xorshiro128p test
./prng/die_xorshiro128p |dieharder -a -g 200