#ifndef ML_H
#define ML_H
#include "all_struct.h"

class roofit_methods{
    //private
    mc_generated vectors;
    estimates est;

    public:
        roofit_methods(mc_generated vec){ vectors = vec ; } ;
        ~roofit_methods(){};
        void ml_unbinned_fit();
        void chi2_binned_fit();

};
#endif //ML_H