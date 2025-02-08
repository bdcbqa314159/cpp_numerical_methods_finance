#include "coding_along.hpp"
#include <cassert>

double risk_neutral_calculator(double U, double D, double R)
{
    bool non_arbitrage = D < R && R < U;
    assert(non_arbitrage && "we must have D < R < U to avoid arbitrage");

    double q = (U - R) / (D - R);
    return q;
}

double h_call(double x, double K)
{
    return (x > K) * (x - K);
}

double h_put(double x, double K)
{
    return (x < K) * (K - x);
}