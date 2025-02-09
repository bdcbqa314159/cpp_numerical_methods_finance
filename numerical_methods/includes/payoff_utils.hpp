#pragma once

double h_call(double x, double K);
double h_put(double x, double K);

double h_digital_call(double x, double K);
double h_digital_put(double x, double K);

double h_double_digital(double x, double K1, double K2);

double h_bull_spread(double x, double K1, double K2);
double h_bear_spread(double x, double K1, double K2);

double h_strangle(double x, double K1, double K2);
double h_butterfly(double x, double K1, double K2);