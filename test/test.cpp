#include <iostream>
#include "numerical_methods"

int main()
{
    double U = 0.9, D = 0.2, R = 0.5;
    double q = risk_neutral_calculator(U, D, R);

    std::cout << "risk neutral q = " << q << "\n";

    D = 0.6;
    q = risk_neutral_calculator(U, D, R);
    std::cout << "risk neutral q = " << q << "\n";

    return 0;
}
