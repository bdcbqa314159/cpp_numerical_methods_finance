#include "payoff_utils.hpp"

double h_call(double x, double K)
{
    return (x > K) * (x - K);
}

double h_put(double x, double K)
{
    return (x < K) * (K - x);
}

double h_digital_call(double x, double K)
{
    return (K < x);
}

double h_digital_put(double x, double K)
{
    return (K > x);
}

double h_double_digital(double x, double K1, double K2)
{
    return (K1 < x && x < K2);
}

double h_bull_spread(double x, double K1, double K2)
{
    if (x <= K1)
        return 0.;

    else if (K1 < x && x < K2)
        return x - K1;

    else
        return K2 - K1;
}

double h_bear_spread(double x, double K1, double K2)
{
    if (x <= K1)
        return K2 - K1;

    else if (K1 < x && x < K2)
        return K2 - x;

    else
        return 0.;
}

double h_strangle(double x, double K1, double K2)
{
    if (x <= K1)
        return K1 - x;
    else if (K1 < x && x < K2)
        return 0.;
    else
        return x - K2;
}

double h_butterfly(double x, double K1, double K2)
{
    double mid_Ks = K1 + 0.5 * (K2 - K1);

    if (K1 < x && x <= mid_Ks)
        return x - K1;
    else if (mid_Ks < x && x < K2)
        return K2 - x;
    else
        return 0.;
}