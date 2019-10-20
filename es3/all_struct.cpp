#ifndef ALL_STRUCT_H
#define ALL_STRUCT_H

struct estimates{ //struct to contain result of fit
    double a_est;
    double b_est;
    double g_est;
};

struct mc_generated{
    std::vector<double> theta;
    std::vector<double> phi;
};

#endif //ALL_STRUCT_H
