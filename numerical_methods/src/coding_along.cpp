#include "coding_along.hpp"
#include <cassert>

double risk_neutral_calculator(double U, double D, double R)
{
    bool non_arbitrage = D < R && R < U;
    assert(non_arbitrage && "we must have D < R < U to avoid arbitrage");

    double q = (R - D) / (U - D);
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

double binomial_price_stock(double S0, double U, double D, int n, int i)
{
    bool price_positive = S0 >= 0.;
    bool percents = (U >= 0. && U <= 1.) && (D >= 0. && D <= 1.);
    bool good_nodes = n >= i && i >= 0;

    assert(price_positive);
    assert(percents);
    assert(good_nodes);

    return S0 * std::pow(1 + U, i) * std::pow(1 + D, n - i);
}