#include <iostream>
#include <iostream>
#include <cmath>
#include <cassert>

using namespace std;

int GetInputData(int &N, std::vector<double> &Ks);

double PriceByCRR(double S0, double U, double D,
                  double R, int N, std::vector<double> &Ks,
                  double (*Payoff)(double z, std::vector<double> &));

double CallPayoff(double z, std::vector<double> &);

double PutPayoff(double z, std::vector<double> &);

double DigitalCallPayoff(double z, std::vector<double> &);

double DigitalPutPayoff(double z, std::vector<double> &);

double DoubleDigitalPayoff(double z, std::vector<double> &);

double RiskNeutProb(double U, double D, double R);

double S(double S0, double U, double D, int n, int i);

int GetInputData(double &S0, double &U, double &D, double &R);

int main()
{
    double S0, U, D, R;

    if (GetInputData(S0, U, D, R) == 1)
        return 1;

    std::vector<double> Ks; // strikes
    int N;                  // steps to expiry

    cout << "Enter call option data:" << endl;
    GetInputData(N, Ks);
    cout << "European call option price = "
         << PriceByCRR(S0, U, D, R, N, Ks, CallPayoff)
         << endl
         << endl;

    cout << "Enter put option data:" << endl;
    GetInputData(N, Ks);
    cout << "European put option price =  "
         << PriceByCRR(S0, U, D, R, N, Ks, PutPayoff)
         << endl
         << endl;

    GetInputData(N, Ks);
    cout << "European digital call option price =  "
         << PriceByCRR(S0, U, D, R, N, Ks, DigitalCallPayoff)
         << endl
         << endl;

    GetInputData(N, Ks);
    cout << "European digital put option price =  "
         << PriceByCRR(S0, U, D, R, N, Ks, DigitalPutPayoff)
         << endl
         << endl;

    return 0;
}

int GetInputData(int &N, std::vector<double> &Ks)
{
    cout << "Enter steps to expiry N: ";
    cin >> N;
    int n{};
    cout << "Enter number of strikes (min 1, max 2):    ";
    cin >> n;
    assert(n == 1 || 2);
    Ks.resize(n);

    if (n == 1)
    {
        cout << "Enter strike price K:    ";
        cin >> Ks[0];
        cout << endl;
    }

    else
    {

        cout << "Enter strike price K1:    ";
        cin >> Ks[0];
        cout << endl;
        cout << "Enter strike price K1:    ";
        cin >> Ks[1];
        cout << endl;
    }

    return 0;
}

double PriceByCRR(double S0, double U, double D,
                  double R, int N, std::vector<double> &Ks,
                  double (*Payoff)(double z, std::vector<double> &))
{
    double q = RiskNeutProb(U, D, R);
    double Price[N + 1];
    for (int i = 0; i <= N; i++)
    {
        Price[i] = Payoff(S(S0, U, D, N, i), Ks);
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

double CallPayoff(double z, std::vector<double> &K)
{
    assert(K.size() == 1);
    if (z > K[0])
        return z - K[0];
    return 0.0;
}

double PutPayoff(double z, std::vector<double> &K)
{
    assert(K.size() == 1);
    if (z < K[0])
        return K[0] - z;
    return 0.0;
}

double DigitalCallPayoff(double z, std::vector<double> &K)
{
    assert(K.size() == 1);
    if (z > K[0])
        return 1.;
    return 0.;
}

double DigitalPutPayoff(double z, std::vector<double> &K)
{
    assert(K.size() == 1);
    if (K[0] > z)
        return 1.;
    return 0.;
}

double DoubleDigitalPayoff(double z, std::vector<double> &K)
{
    assert(K.size() == 2 && K[0] < K[1]);
    if (z < K[1] && z > K[0])
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
