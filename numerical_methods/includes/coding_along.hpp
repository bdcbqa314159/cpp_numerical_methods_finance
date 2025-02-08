#pragma once
#include <iostream>

// putting my coded notes here the structure will organised later

double risk_neutral_calculator(double U, double D, double R);
double h_call(double x, double K);
double h_put(double x, double K);

void binomial_parameters_checker(const double &S0, const double &U, const double &D, const int &n, const int &i);
double binomial_price_stock(double S0, double U, double D, int n, int i);
