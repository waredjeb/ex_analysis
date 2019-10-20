#include "TRandom3.h"
#include "prng.h"
#include <climits>
#include <iostream>



int main(int argc, char** argv){

    TRandom3* r = new TRandom3();

    randtype rnd_num = r->Uniform(-UINT_MAX,UINT_MAX);

    std::cout.write(reinterpret_cast<char*>(&rnd_num), sizeof rnd_num);

    while(1){

    rnd_num = r->Uniform(-UINT_MAX,UINT_MAX);
    //std::cout << rnd_num << std::endl;

    std::cout.write(reinterpret_cast<char*>(&rnd_num), sizeof rnd_num);

    }


    return 0;
}
