#include <iostream>
#include <iomanip>
using namespace std;

template <typename Type>
class BinLattice
{
private:
    int N;
    vector<vector<Type>> Lattice;

public:
    void SetN(int N_)
    {
        N = N_;
        Lattice.resize(N + 1);
        for (int n = 0; n <= N; n++)
            Lattice[n].resize(n + 1);
    }
    void SetNode(int n, int i, Type x)
    {
        Lattice[n][i] = x;
    }
    Type GetNode(int n, int i)
    {
        return Lattice[n][i];
    }
    void Display()
    {
        cout << setiosflags(ios::fixed)
             << setprecision(3);
        for (int n = 0; n <= N; n++)
        {
            for (int i = 0; i <= n; i++)
                cout << setw(7) << GetNode(n, i);
            cout << endl;
        }
        cout << endl;
    }
};

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

class Option
{
private:
    int N; // steps to expiry

public:
    void SetN(int N_) { N = N_; }
    int GetN() { return N; }
    virtual double Payoff(double z) = 0;
};

class EurOption : public virtual Option
{
public:
    // pricing European option
    // double PriceByCRR(BinModel Model);
};

class AmOption : public virtual Option
{
public:
    // pricing American option
    // double PriceBySnell(BinModel Model,
    //                     BinLattice<double> &PriceTree,
    //                     BinLattice<bool> &StoppingTree);
};

class Call : public EurOption, public AmOption
{
private:
    double K; // strike price

public:
    void SetK(double K_) { K = K_; }
    int GetInputData();
    double Payoff(double z);
};

class Put : public EurOption, public AmOption
{
private:
    double K; // strike price

public:
    void SetK(double K_) { K = K_; }
    int GetInputData();
    double Payoff(double z);
};

template <typename Option>
double PriceByCRR(BinModel Model, Option theOption)
{
    double q = Model.RiskNeutProb();
    int N = theOption.GetN();
    vector<double> Price(N + 1);
    for (int i = 0; i <= N; i++)
    {
        Price[i] = theOption.Payoff(Model.S(N, i));
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

template <typename Option>
double PriceBySnell(BinModel Model,
                    BinLattice<double> &PriceTree,
                    BinLattice<bool> &StoppingTree, Option theOption)
{
    double q = Model.RiskNeutProb();
    int N = theOption.GetN();
    PriceTree.SetN(N);
    StoppingTree.SetN(N);
    double ContVal;
    for (int i = 0; i <= N; i++)
    {
        PriceTree.SetNode(N, i, theOption.Payoff(Model.S(N, i)));
        StoppingTree.SetNode(N, i, 1);
    }
    for (int n = N - 1; n >= 0; n--)
    {
        for (int i = 0; i <= n; i++)
        {
            ContVal = (q * PriceTree.GetNode(n + 1, i + 1) + (1 - q) * PriceTree.GetNode(n + 1, i)) / (1 + Model.GetR());
            PriceTree.SetNode(n, i, theOption.Payoff(Model.S(n, i)));
            StoppingTree.SetNode(n, i, 1);
            if (ContVal > PriceTree.GetNode(n, i))
            {
                PriceTree.SetNode(n, i, ContVal);
                StoppingTree.SetNode(n, i, 0);
            }
            else if (PriceTree.GetNode(n, i) == 0.0)
            {
                StoppingTree.SetNode(n, i, 0);
            }
        }
    }
    return PriceTree.GetNode(0, 0);
}

int main()
{
    BinModel Model;

    if (Model.GetInputData() == 1)
        return 1;

    Put Option;
    Option.GetInputData();
    BinLattice<double> PriceTree;
    BinLattice<bool> StoppingTree;
    PriceBySnell(Model, PriceTree, StoppingTree, Option);
    cout << "American put prices:" << endl
         << endl;
    PriceTree.Display();
    cout << "American put exercise policy:"
         << endl
         << endl;
    StoppingTree.Display();
    return 0;
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
    return 0;
}

double Put::Payoff(double z)
{
    if (z < K)
        return K - z;
    return 0.0;
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
