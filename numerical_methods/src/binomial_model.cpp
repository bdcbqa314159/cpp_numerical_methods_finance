#include "binomial_model.hpp"

void BinomialModel::setS0(double S0)
{
    initial_stock_price_checker(S0);
    this->S0 = S0;
}

void BinomialModel::setDirections(double U, double D)
{
    direction_checker(U, D);
    this->U = U;
    this->D = D;
}

void BinomialModel::setR(double R)
{
    risk_neutral_checker(U, D, R);
    this->R = R;
}

void BinomialModel::setParameters(double S0, double U, double D, double R)
{
    setS0(S0);
    setDirections(U, D);
    setR(R);
}

double BinomialModel::riskNeutralProbability() const
{
    return risk_neutral_calculator(U, D, R);
}

double BinomialModel::S(int n, int i) const
{
    return binomial_price_stock(S0, U, D, n, i);
}