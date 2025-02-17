#include <iostream>

#include <vector>
#include <cstdlib>
#include <ctime>

typedef std::vector<double> SamplePath;
const double pi = 4.0 * atan(1.0);

class BSModel
{
public:
    double S0, r, sigma;
    BSModel(double S0_, double r_, double sigma_)
    {
        S0 = S0_;
        r = r_;
        sigma = sigma_;
        srand(time(NULL));
    }
    void GenerateSamplePath(double T, int m, SamplePath &S);
};

class Option
{
public:
    double T, zl, zu;
    virtual double Payoff(double z) = 0;
    virtual double UpperBdCond(BSModel *PtrModel, double t) = 0;
    virtual double LowerBdCond(BSModel *PtrModel, double t) = 0;
};

class Put : public Option
{
public:
    double K;
    Put(double K_, double T_, double zl_, double zu_)
    {
        K = K_;
        T = T_;
        zl = zl_;
        zu = zu_;
    }
    double Payoff(double z);
    double UpperBdCond(BSModel *PtrModel, double t);
    double LowerBdCond(BSModel *PtrModel, double t);
};

int main()
{
    std::cout << "Extending the code to have calls... \n";
    return 0;
}

double Gauss()
{
    double U1 = (rand() + 1.0) / (RAND_MAX + 1.0);
    double U2 = (rand() + 1.0) / (RAND_MAX + 1.0);
    return sqrt(-2.0 * log(U1)) * cos(2.0 * pi * U2);
}

void BSModel::GenerateSamplePath(double T, int m, SamplePath &S)
{
    double St = S0;
    for (int k = 0; k < m; k++)
    {
        S[k] = St * exp((r - sigma * sigma * 0.5) * (T / m) + sigma * sqrt(T / m) * Gauss());
        St = S[k];
    }
}

double Put::Payoff(double z)
{
    if (K < z)
        return 0.0;
    return K - z;
}

double Put::UpperBdCond(BSModel *PtrModel, double t)
{
    return 0.0;
}

double Put::LowerBdCond(BSModel *PtrModel, double t)
{
    return K * exp(-PtrModel->r * (T - t));
}