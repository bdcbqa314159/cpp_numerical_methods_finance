#ifndef PathDepOption05_h
#define PathDepOption05_h

#include "bs_model02.hpp"

class PathDepOption
{
public:
    double T;
    int m;
    double PriceByMC(BSModel Model, long N);
    virtual double Payoff(SamplePath &S) = 0;
};

class ArthmAsianCall : public PathDepOption
{
public:
    double K;
    ArthmAsianCall(double T_, double K_, int m_)
    {
        T = T_;
        K = K_;
        m = m_;
    }
    double Payoff(SamplePath &S);
};

#endif