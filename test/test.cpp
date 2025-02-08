#include <iostream>
#include "numerical_methods"

int main()
{
    double U = 0.9, D = 0.2, R = 0.5, q{};
    q = risk_neutral_calculator(U, D, R);

    std::cout << "risk neutral q = " << q << "\n";

    // D = 0.6;
    // q = risk_neutral_calculator(U, D, R);
    // std::cout << "risk neutral q = " << q << "\n";

    int n = 3, i = 2;
    D = 0.1;
    double S0 = 100.;
    double S = binomial_price_stock(S0, U, D, n, i);

    std::cout << "binomial price stock S = " << S << "\n";
    return 0;
}
