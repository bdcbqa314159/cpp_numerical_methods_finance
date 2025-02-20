#ifndef BSEqLCP_h
#define BSEqLCP_h

#include "lcp.hpp"
#include "bs_model01.hpp"
#include "option.hpp"
#include "bs_eq.hpp"

class BSEqLCP : public LCP, public BSEq
{
public:
    BSEqLCP(BSModel *PtrModel, Option *PtrOption)
        : BSEq(PtrModel, PtrOption) { PtrPDE = this; }

    double g(double t, double z)
    {
        return PtrOption->Payoff(z);
    }
};

#endif
