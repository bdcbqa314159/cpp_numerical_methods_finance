/*
Include checking for input data integrity in GetInputData() function in Options01.cpp.
You want to ensure that 0<K and 0<N.
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

long int my_factorial(int n)
{
    if (n == 0 || n == 1)
        return 1;
    long result = 1;
    for (int j = 1; j <= n; ++j)
        result *= j;
    return result;
}

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

    // checking for arbitrage
    if (K <= 0. || N <= 0)
    {
        cout << "Bad maturity or strike" << endl;
        cout << "Terminating program" << endl;
        return 1;
    }
    return 0;
}

double PriceByCRR(double S0, double U, double D,
                  double R, int N, double K)
{

    double q = RiskNeutProb(U, D, R);

    double discount = 1. / (std::pow(1 + R, N));
    long int N_factorial = my_factorial(N);
    double p = 1 - q;
    double result = 0.;

    for (int i = 0; i <= N; i++)
    {
        long int i_factorial = my_factorial(i);
        long int N_i_factorial = my_factorial(N - i);

        double q_i = std::pow(q, i);
        double p_N_i = std::pow(p, N - i);

        double call_S_N_i = CallPayoff(S(S0, U, D, N, i), K);

        double numerator = N_factorial * q_i * p_N_i * call_S_N_i;
        double denominator = i_factorial * N_i_factorial;
        result += (numerator / denominator);
    }
    result *= discount;
    return result;
}

double CallPayoff(double z, double K)
{
    if (z > K)
        return z - K;
    return 0.0;
}