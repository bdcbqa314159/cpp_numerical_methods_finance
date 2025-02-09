#include "binomial_model.hpp"

void BinomialModel::setParameters(double S0, double U, double D, double R)
{
    risk_neutral_checker(U, D, R);
    this->S0 = S0;
    this->U = U;
    this->D = D;
    this->R = R;
}

double BinomialModel::riskNeutralProbability() const
{
    return risk_neutral_calculator(U, D, R);
}

double BinomialModel::S(int n, int i) const
{
    return binomial_price_stock(S0, U, D, n, i);
}