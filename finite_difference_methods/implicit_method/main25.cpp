#include <iostream>
#include "bs_model01.hpp"
#include "option.hpp"
#include "bs_eq.hpp"
#include "cn_method.hpp"

int main()
{
    double S0 = 100.0, r = 0.05, sigma = 0.2;
    BSModel Model(S0, r, sigma);

    double T = 1. / 12., K = 100.0, zl = 0.0, zu = 2.0 * S0;
    Put EuropeanPut(K, T, zl, zu);

    int imax = 200, jmax = 2000;

    BSEq BSPDE(&Model, &EuropeanPut);

    CNMethod Method(&BSPDE, imax, jmax);
    Method.SolvePDE();
    cout << "Price = " << Method.v(0.0, S0) << endl;

    return 0;
}
