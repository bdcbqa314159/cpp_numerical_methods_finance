#include "coding_along.hpp"

double EurOption::price_by_CRR(BinomialModel model)
{
    double q = model.riskNeutralProbability();
    double p = 1 - q;
    int N = get_N();
    std::vector<double> price(N + 1);
    double discount = 1. / (1 + model.getR());

    for (int i = 0; i <= N; ++i)
        price[i] = payoff(model.S(N, i));

    for (int n = N - 1; n >= 0; --n)
    {
        for (int i = 0; i <= n; ++i)
            price[i] = discount * (q * price[i + 1] + p * price[i]);
    }
    return price[0];
}

double Call::payoff(double x) const
{
    return (x > K) * (x - K);
}

double Put::payoff(double x) const
{
    return (x < K) * (K - x);
}

double DigitalCall::payoff(double x) const
{
    return (x > K) * 1.;
}

double DigitalPut::payoff(double x) const
{
    return (x < K) * 1.;
}

double BullSpread::payoff(double x) const
{
    if (x <= K1)
        return 0.;
    else if (K1 < x && x < K2)
    {
        return x - K1;
    }
    else
    {
        return K2 - K1;
    }
}

double BearSpread::payoff(double x) const
{
    if (K2 < x)
        return 0.;
    else if (K1 < x && x < K2)
    {
        return K2 - x;
    }
    else
    {
        return K2 - K1;
    }
}

double DoubleDigital::payoff(double x) const
{
    if (K1 < x && x < K2)
        return 1.;
    return 0.;
}

double Strangle::payoff(double x) const
{
    if (x <= K1)
        return K1 - x;
    else if (K1 < x && x <= K2)
    {
        return 0.;
    }
    else
    {
        return x - K2;
    }
}

double Butterfly::payoff(double x) const
{
    double midK = K1 + 0.5 * (K2 - K1);

    if (K1 < x && x <= midK)
        return x - K1;

    else if (midK < x && x < K2)
    {
        return K2 - x;
    }
    else
    {
        return 0.;
    }
}