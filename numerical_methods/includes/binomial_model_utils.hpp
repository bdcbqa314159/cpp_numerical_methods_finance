#pragma once

void binomial_parameters_checker(const double &S0, const double &U, const double &D, const int &n, const int &i);

double binomial_price_stock(double S0, double U, double D, int n, int i);

void risk_neutral_checker(double &U, double &D, double &R);

double risk_neutral_calculator(double U, double D, double R);