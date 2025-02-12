#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace std;

const double pi = 4.0 * atan(1.0);
typedef vector<double> SamplePath;

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

class PathDepOption
{
public:
    double T;
    int m;
    double PriceByMC(BSModel Model, long N);
    virtual double Payoff(SamplePath &S) = 0;
};

class ArthmAsianCall : public PathDepOption
{
public:
    double K;
    ArthmAsianCall(double T_, double K_, int m_)
    {
        T = T_;
        K = K_;
        m = m_;
    }
    double Payoff(SamplePath &S);
};

class EurCall : public PathDepOption
{
public:
    double K;
    EurCall(double T_, double K_)
    {
        T = T_;
        K = K_;
        m = 1;
    }
    double Payoff(SamplePath &S);
};

class EurPut : public PathDepOption
{
public:
    double K;
    EurPut(double T_, double K_)
    {
        T = T_;
        K = K_;
        m = 1;
    }
    double Payoff(SamplePath &S);
};

int main()
{
    double S0 = 100.0, r = 0.03, sigma = 0.2;
    BSModel Model(S0, r, sigma);

    double T = 1.0 / 12.0, K = 100.0;
    int m = 30;
    ArthmAsianCall Option(T, K, m);

    long N = 30000;
    cout << "Asian Call Price = "
         << Option.PriceByMC(Model, N) << endl;

    EurCall CallOption(T, K);
    EurPut PutOption(T, K);

    cout << "Eur Call Price = "
         << CallOption.PriceByMC(Model, N) << endl;

    cout << "Eur Put Price = "
         << PutOption.PriceByMC(Model, N) << endl;

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

double PathDepOption::PriceByMC(BSModel Model, long N)
{
    double H = 0.0;
    SamplePath S(m);
    for (long i = 0; i < N; i++)
    {
        Model.GenerateSamplePath(T, m, S);
        H = (i * H + Payoff(S)) / (i + 1.0);
    }
    return exp(-Model.r * T) * H;
}

double ArthmAsianCall::Payoff(SamplePath &S)
{
    double Ave = 0.0;
    for (int k = 0; k < m; k++)
        Ave = (k * Ave + S[k]) / (k + 1.0);
    if (Ave < K)
        return 0.0;
    return Ave - K;
}

double EurCall::Payoff(SamplePath &S)
{
    if (S[0] > K)
        return S[0] - K;
    return 0.;
}

double EurPut::Payoff(SamplePath &S)
{
    if (S[0] < K)
        return K - S[0];
    return 0.;
}