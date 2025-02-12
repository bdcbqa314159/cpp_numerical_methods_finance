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
    double T, Price, PricingError, delta, gamma;
    int m;
    double PriceByMC(BSModel Model, long N, double epsilon);
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

void Rescale(SamplePath &S, double x);

int main()
{
    double S0 = 100.0, r = 0.03, sigma = 0.2, epsilon = 0.0001;
    BSModel Model(S0, r, sigma);

    double T = 1.0 / 12.0, K = 100.0;
    int m = 30;
    ArthmAsianCall Option(T, K, m);

    long N = 30000;
    cout << "Asian Call Price = "
         << Option.PriceByMC(Model, N, epsilon) << endl
         << "Pricing Error = " << Option.PricingError << endl
         << "delta = " << Option.delta << endl
         << "gamma = " << Option.gamma << endl;

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

void Rescale(SamplePath &S, double x)
{
    int m = S.size();
    for (int i = 0; i < m; ++i)
        S[i] *= x;
}

double PathDepOption::PriceByMC(BSModel Model, long N, double epsilon)
{
    double H = 0.0, Hsq = 0.0, Heps = 0.0, Hmeps = 0.0;
    SamplePath S(m);
    for (long i = 0; i < N; i++)
    {
        Model.GenerateSamplePath(T, m, S);
        H = (i * H + Payoff(S)) / (i + 1.0);
        Hsq = (i * Hsq + pow(Payoff(S), 2)) / (i + 1.);
        Rescale(S, 1. + epsilon);
        Heps = (i * Heps + Payoff(S)) / (i + 1.);
        Rescale(S, (1. - epsilon) / (1. + epsilon));
        Hmeps = (i * Hmeps + Payoff(S)) / (i + 1.);
    }
    Price = exp(-Model.r * T) * H;
    PricingError = exp(-Model.r * T) * sqrt(Hsq - H * H) / sqrt(N - 1.);
    delta = exp(-Model.r * T) * (Heps - H) / (Model.S0 * epsilon);
    gamma = exp(-Model.r * T) * (Heps - 2 * H + Hmeps) / ((Model.S0 * epsilon) * (Model.S0 * epsilon));
    return Price;
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