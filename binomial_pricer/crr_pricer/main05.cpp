#include "bin_model01.hpp"
#include "options01.hpp"
#include <iostream>
#include <cmath>
using namespace std;

int main()
{
    double S0, U, D, R;

    if (GetInputData(S0, U, D, R) == 1)
        return 1;

    double K; // strike price
    int N;    // steps to expiry

    cout << "Enter call option data:" << endl;
    GetInputData(N, K);
    cout << "European call option price = "
         << PriceByCRR(S0, U, D, R, N, K)
         << endl
         << endl;

    return 0;
}
