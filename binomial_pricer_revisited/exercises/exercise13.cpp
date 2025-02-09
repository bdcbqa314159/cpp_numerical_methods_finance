#include <iostream>
#include <cmath>
#include <cassert>
using namespace std;

class BinModel
{
private:
    double S0;
    double U;
    double D;
    double R;

public:
    // computing risk-neutral probability
    double RiskNeutProb();

    // computing the stock price at node n,i
    double S(int n, int i);

    // inputting, displaying and checking model data
    int GetInputData();

    double GetR();
};

class EurOption
{
private:
    int N; // steps to expiry
public:
    void SetN(int N_) { N = N_; }
    // Payoff defined to return 0.0
    // for pedagogical purposes.
    // To use a pure virtual function replace by
    // virtual double Payoff(double z)=0;
    virtual double Payoff(double z) { return 0.0; }
    // pricing European option
    double PriceByCRR(BinModel Model);
};

class Call : public EurOption
{
private:
    double K; // strike price
public:
    void SetK(double K_) { K = K_; }
    int GetInputData();
    double Payoff(double z);
};

class Put : public EurOption
{
private:
    double K; // strike price
public:
    void SetK(double K_) { K = K_; }
    int GetInputData();
    double Payoff(double z);
};

class DoubleStrikeOption : public EurOption
{
private:
    double K1, K2; // strike price
public:
    void SetK(double K1, double K2)
    {
        this->K1 = K1;
        this->K2 = K2;
    }
    int GetInputData();
    double Payoff(double z) { return 0.; }
};

class BullSpread : public DoubleStrikeOption
{
private:
    double K1, K2; // strike price
public:
    double Payoff(double z);
};

class BearSpread : public DoubleStrikeOption
{
private:
    double K1, K2; // strike price
public:
    double Payoff(double z);
};

int main()
{
    BinModel Model;

    if (Model.GetInputData() == 1)
        return 1;

    Call Option1;
    Option1.GetInputData();
    cout << "European call option price = "
         << Option1.PriceByCRR(Model)
         << endl
         << endl;

    Put Option2;
    Option2.GetInputData();
    cout << "European put option price = "
         << Option2.PriceByCRR(Model)
         << endl
         << endl;

    BearSpread Option3;
    Option3.GetInputData();
    cout << "European bear spread option price = "
         << Option3.PriceByCRR(Model)
         << endl
         << endl;

    BullSpread Option4;
    Option4.GetInputData();
    cout << "European bull spread option price = "
         << Option4.PriceByCRR(Model)
         << endl
         << endl;

    return 0;
}

double EurOption::PriceByCRR(BinModel Model)
{
    double q = Model.RiskNeutProb();
    double Price[N + 1];
    for (int i = 0; i <= N; i++)
    {
        Price[i] = Payoff(Model.S(N, i));
    }
    for (int n = N - 1; n >= 0; n--)
    {
        for (int i = 0; i <= n; i++)
        {
            Price[i] = (q * Price[i + 1] + (1 - q) * Price[i]) / (1 + Model.GetR());
        }
    }
    return Price[0];
}

int Call::GetInputData()
{
    cout << "Enter call option data:" << endl;
    int N;
    cout << "Enter steps to expiry N: ";
    cin >> N;
    SetN(N);
    cout << "Enter strike price K:    ";
    cin >> K;
    cout << endl;
    SetK(K);
    return 0;
}

double Call::Payoff(double z)
{
    if (z > K)
        return z - K;
    return 0.0;
}

int Put::GetInputData()
{
    cout << "Enter put option data:" << endl;
    int N;
    cout << "Enter steps to expiry N: ";
    cin >> N;
    SetN(N);
    cout << "Enter strike price K:    ";
    cin >> K;
    cout << endl;
    SetK(K);
    return 0;
}

double Put::Payoff(double z)
{
    if (z < K)
        return K - z;
    return 0.0;
}

int DoubleStrikeOption::GetInputData()
{
    cout << "Enter double strike option data:" << endl;
    int N;
    cout << "Enter steps to expiry N: ";
    cin >> N;
    SetN(N);
    cout << "Enter strike price K1:    ";
    cin >> K1;
    cout << "Enter strike price K2:    ";
    cin >> K2;
    cout << endl;

    assert(K1 < K2);
    SetK(K1, K2);
    return 0;
}

double BullSpread::Payoff(double z)
{
    if (z <= K1)
        return 0.;
    else if (K1 < z && z < K2)
    {
        return z - K1;
    }
    else
    {
        return K2 - K1;
    }
}

double BearSpread::Payoff(double z)
{
    if (K2 < z)
        return 0.;
    else if (K1 < z && z < K2)
    {
        return K2 - z;
    }
    else
    {
        return K2 - K1;
    }
}

double BinModel::RiskNeutProb()
{
    return (R - D) / (U - D);
}

double BinModel::S(int n, int i)
{
    return S0 * pow(1 + U, i) * pow(1 + D, n - i);
}

int BinModel::GetInputData()
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

double BinModel::GetR()
{
    return R;
}
