#include <iostream>
#include <iostream>
#include <cmath>
using namespace std;

int GetInputData(int &N, double &K);

double PriceByCRR(double S0, double U, double D, double R, int N, double K, double (*Payoff)(double z, double K));

double CallPayoff(double z, double K);

double PutPayoff(double z, double K);

double DigitalPutPayoff(double z, double K);

double RiskNeutProb(double U, double D, double R);

double S(double S0, double U, double D, int n, int i);

int GetInputData(double &S0, double &U, double &D, double &R);

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
         << PriceByCRR(S0, U, D, R, N, K, CallPayoff)
         << endl
         << endl;

    cout << "Enter put option data:" << endl;
    GetInputData(N, K);
    cout << "European put option price =  "
         << PriceByCRR(S0, U, D, R, N, K, PutPayoff)
         << endl
         << endl;

    GetInputData(N, K);
    cout << "European digital put option price =  "
         << PriceByCRR(S0, U, D, R, N, K, DigitalPutPayoff)
         << endl
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
                  double R, int N, double K,
                  double (*Payoff)(double z, double K))
{
    double q = RiskNeutProb(U, D, R);
    double Price[N + 1];
    for (int i = 0; i <= N; i++)
    {
        Price[i] = Payoff(S(S0, U, D, N, i), K);
    }
    for (int n = N - 1; n >= 0; n--)
    {
        for (int i = 0; i <= n; i++)
        {
            Price[i] = (q * Price[i + 1] + (1 - q) * Price[i]) / (1 + R);
        }
    }
    return Price[0];
}

double CallPayoff(double z, double K)
{
    if (z > K)
        return z - K;
    return 0.0;
}

double PutPayoff(double z, double K)
{
    if (z < K)
        return K - z;
    return 0.0;
}

double DigitalPutPayoff(double z, double K)
{
    if (K > z)
        return 1.;
    return 0.;
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
