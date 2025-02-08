/*
Rewrite the loops in PriceByCRR() in Options01.cpp using while instead
of for. You may need to read about the while loop from a c++ manual.
*/

/*
To compare:

Enter S0: 100
Enter U:  0.9
Enter D:  0.2
Enter R:  0.5

Input data checked
There is no arbitrage

Enter call option data:
Enter steps to expiry N: 10
Enter strike price K:    30

European call option price = 99.4798
*/

#include <iostream>
#include <cmath>
using namespace std;

double RiskNeutProb(double U, double D, double R);

double S(double S0, double U, double D, int n, int i);

int GetInputData(double &S0, double &U, double &D, double &R);

int GetInputData(int &N, double &K);

double PriceByCRR(double S0, double U, double D, double R, int N, double K);

double CallPayoff(double z, double K);

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

double RiskNeutProb(double U, double D, double R)
{
    return (R - D) / (U - D);
}

double S(double S0, double U, double D, int n, int i)
{
    return S0 * pow(1 + U, i) * pow(1 + D, n - i);
}

int GetInputData(double &S0,
                 double &U, double &D, double &R)
{
    // entering data
    cout << "Enter S0: ";
    cin >> S0;
    cout << "Enter U:  ";
    cin >> U;
    cout << "Enter D:  ";
    cin >> D;
    cout << "Enter R:  ";
    cin >> R;
    cout << endl;

    // making sure that 0<S0, -1<D<U, -1<R
    if (S0 <= 0.0 || U <= -1.0 || D <= -1.0 || U <= D || R <= -1.0)
    {
        cout << "Illegal data ranges" << endl;
        cout << "Terminating program" << endl;
        return 1;
    }

    // checking for arbitrage
    if (R >= U || R <= D)
    {
        cout << "Arbitrage exists" << endl;
        cout << "Terminating program" << endl;
        return 1;
    }

    cout << "Input data checked" << endl;
    cout << "There is no arbitrage" << endl
         << endl;

    return 0;
}

int GetInputData(int &N, double &K)
{
    cout << "Enter steps to expiry N: ";
    cin >> N;
    cout << "Enter strike price K:    ";
    cin >> K;
    cout << endl;

    return 0;
}

double PriceByCRR(double S0, double U, double D,
                  double R, int N, double K)
{
    double q = RiskNeutProb(U, D, R);
    double Price[N + 1];

    int i = 0, n = N - 1;

    while (i <= N)
    {
        Price[i] = CallPayoff(S(S0, U, D, N, i), K);
        ++i;
    }

    while (n >= 0)
    {
        i = 0;
        while (i <= n)
        {
            Price[i] = (q * Price[i + 1] + (1 - q) * Price[i]) / (1 + R);
            ++i;
        }

        --n;
    }

    return Price[0];
}

double CallPayoff(double z, double K)
{
    if (z > K)
        return z - K;
    return 0.0;
}