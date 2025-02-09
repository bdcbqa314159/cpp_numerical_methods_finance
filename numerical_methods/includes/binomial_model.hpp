#pragma once
#include "binomial_model_utils.hpp"

class BinomialModel
{
private:
    double S0 = 100., U = 0.9, D = 0.1, R = 0.5;

public:
    BinomialModel() = default;
    BinomialModel(double S0, double U, double D, double R) : S0(S0), U(U), D(D), R(R)
    {
        risk_neutral_checker(U, D, R);
    }

    double getR() const { return R; };
    void setParameters(double S0, double U, double D, double R);
    double riskNeutralProbability() const;
    double S(int n, int i) const;
};