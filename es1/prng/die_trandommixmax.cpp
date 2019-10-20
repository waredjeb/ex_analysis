#include "TRandomMixMax.h"
#include "prng.h"
#include <climits>
#include <iostream>



int main(int argc, char** argv){

    TRandomMixMax* r = new TRandomMixMax();

    randtype rnd_num = r->Uniform(-ULONG_MAX,ULONG_MAX);

    std::cout.write(reinterpret_cast<char*>(&rnd_num), sizeof rnd_num);

    while(1){

    rnd_num = r->Uniform(-ULONG_MAX,ULONG_MAX);
    //std::cout << rnd_num << std::endl;

    std::cout.write(reinterpret_cast<char*>(&rnd_num), sizeof rnd_num);

    }


    return 0;
}