#ifndef LCP_h
#define LCP_h

#include "parab_pde.hpp"

class LCP
{
public:
    ParabPDE *PtrPDE;
    virtual double g(double t, double x) = 0;
};

#endif
