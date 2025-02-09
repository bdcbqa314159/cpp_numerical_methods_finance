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