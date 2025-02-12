#ifndef PathDepOption01_h
#define PathDepOption01_h

#include "bs_model01.hpp"

class PathDepOption
{
public:
    double T, Price, PricingError, delta;
    int m;
    double PriceByMC(BSModel Model, long N, double epsilon);
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

void Rescale(SamplePath &S, double x);

#endif
