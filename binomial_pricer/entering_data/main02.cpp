#include <iostream>

int main()
{
    double S0, U, D, R;

    std::cout << "Enter S0 ";
    std::cin >> S0;
    std::cout << "Enter U ";
    std::cin >> U;
    std::cout << "Enter D ";
    std::cin >> D;
    std::cout << "Enter R ";
    std::cin >> R;
    std::cout << "\n";

    // making sure that 0<S0, -1<D<U, -1<R
    if (S0 <= 0.0 || U <= -1.0 || D <= -1.0 || U <= D || R <= -1.0)
    {
        std::cout << "Illegal data ranges" << std::endl;
        std::cout << "Terminating program" << std::endl;
        return 1;
    }

    // checking for arbitrage
    if (R >= U || R <= D)
    {
        std::cout << "Arbitrage exists" << std::endl;
        std::cout << "Terminating program" << std::endl;
        return 1;
    }

    std::cout << "Input data checked" << std::endl;
    std::cout << "There is no arbitrage" << std::endl
              << std::endl;
    // compute risk-neutral probability
    std::cout << " q=" << (R - D) / (U - D) << std::endl;
    // compute stock price at node n=3,i=2
    int n = 3;
    int i = 2;

    std::cout << " n=" << n << std::endl;
    std::cout << " i=" << i << std::endl;
    std::cout << "S(n,i )=" << S0 * pow(1 + U, i) * pow(1 + D, n - i) << std::endl;

    return 0;
}