#include "coding_along.hpp"
#include "payoff_utils.hpp"
#include <cassert>

void risk_neutral_checker(double &U, double &D, double &R)
{
    bool non_arbitrage = D < R && R < U;
    assert(non_arbitrage);
}

double risk_neutral_calculator(double U, double D, double R)
{
    risk_neutral_checker(U, D, R);

    double q = (R - D) / (U - D);
    return q;
}

void binomial_parameters_checker(const double &S0, const double &U, const double &D, const int &n, const int &i)
{
    bool price_positive = S0 >= 0.;
    bool percents = (U >= 0. && U <= 1.) && (D >= 0. && D <= 1.);
    bool good_nodes = n >= i && i >= 0;

    assert(price_positive);
    assert(percents);
    assert(good_nodes);
}

double binomial_price_stock(double S0, double U, double D, int n, int i)
{
    binomial_parameters_checker(S0, U, D, n, i);
    double S_n_i{};
    S_n_i = S0 * std::pow(1 + U, i) * std::pow(1 + D, n - i);
    return S_n_i;
}