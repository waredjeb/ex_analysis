#include "prng.h"
#include <climits>
#include <iostream>


int main(int argc, char** argv){

    prng gen;

    while(1){

    randtype rnd_num = gen.randInt();
    std::cout.write(reinterpret_cast<char*>(&rnd_num), sizeof rnd_num);

    }
    
    return 0;
}