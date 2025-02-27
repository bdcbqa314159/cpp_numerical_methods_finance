#include <iostream>
#include "bs_model01.hpp"
#include "option.hpp"
#include "bs_eq.hpp"
#include "explicit_method.hpp"

int main()
{
    double S0 = 100.0, r = 0.05, sigma = 0.2;
    BSModel Model(S0, r, sigma);

    double K = 100.0, T = 1. / 12., zl = 0.0, zu = 2.0 * S0;
    Put EuropeanPut(K, T, zl, zu);

    BSEq BSPDE(&Model, &EuropeanPut);

    int imax = 3000, jmax = 1000;
    ExplicitMethod Method(&BSPDE, imax, jmax);

    Method.SolvePDE();

    cout << "Price = " << Method.v(0.0, S0) << endl;

    return 0;
}
