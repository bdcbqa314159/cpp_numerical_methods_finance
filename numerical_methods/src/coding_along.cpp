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

double AmOption::price_by_Snell(BinomialModel model, BinomialLattice<double> &priceTree, BinomialLattice<bool> &stoppingTree)
{
    double q = model.riskNeutralProbability();
    double p = 1 - q;
    int N = get_N();
    priceTree.set_N(N);
    stoppingTree.set_N(N);
    double discount = 1. / (1 + model.getR());

    for (int i = 0; i <= N; ++i)
    {
        priceTree.set_node(N, i, payoff(model.S(N, i)));
        stoppingTree.set_node(N, i, true);
    }

    for (int n = N - 1; n >= 0; --n)
    {
        for (int i = 0; i <= n; ++i)
        {
            double H1 = priceTree.get_node(n + 1, i);
            double H2 = priceTree.get_node(n + 1, i + 1);
            double cont_val = discount * (p * H1 + q * H2);

            priceTree.set_node(n, i, payoff(model.S(n, i)));
            stoppingTree.set_node(n, i, true);

            double H0 = priceTree.get_node(n, i);

            if (cont_val > H0)
            {
                priceTree.set_node(n, i, cont_val);
                stoppingTree.set_node(n, i, 0);
            }

            else if (H0 < std::numeric_limits<double>::epsilon())
                stoppingTree.set_node(n, i, false);
        }
    }
    return priceTree.get_node(0, 0);
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