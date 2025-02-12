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
    double T, Price, PricingError, delta, rho, vega, theta, gamma;
    int m;
    virtual double Payoff(SamplePath &S) = 0;
    double PriceByMC(BSModel Model, long N,
                     double epsilon);
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

void GetZ(vector<double> &Z, SamplePath &S, double S0,
          double r, double sigma, double dt);

void Rescale(SamplePath &S, vector<double> &Z, double S0,
             double r, double sigma, double dt);

int main()
{
    double S0 = 100., r = 0.03, sigma = 0.2;
    BSModel Model(S0, r, sigma);

    double T = 0.5, K = 100.0;
    int m = 30;

    long N = 30000;
    double epsilon = 0.001;
    ArthmAsianCall Option1(T, K, m);
    Option1.PriceByMC(Model, N, epsilon);
    cout << "Arithmetic Call:" << endl;
    cout << " Option Price = " << Option1.Price << endl
         << "Pricing Error = " << Option1.PricingError << endl
         << "        delta = " << Option1.delta << endl
         << "          rho = " << Option1.rho << endl
         << "         vega = " << Option1.vega << endl
         << "        theta = " << Option1.theta << endl
         << "        gamma = " << Option1.gamma << endl;

    EurCall Option2(T, K);
    cout << "European Call:" << endl;
    Option2.PriceByMC(Model, N, epsilon);
    cout << " Option Price = " << Option2.Price << endl
         << "Pricing Error = " << Option2.PricingError << endl
         << "        delta = " << Option2.delta << endl
         << "          rho = " << Option2.rho << endl
         << "         vega = " << Option2.vega << endl
         << "        theta = " << Option2.theta << endl
         << "        gamma = " << Option2.gamma << endl;

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

void GetZ(vector<double> &Z, SamplePath &S, double S0,
          double r, double sigma, double dt)
{
    int m = S.size();
    Z[0] = (log(S[0] / S0) - (r - sigma * sigma / 2.0) * dt) / (sigma * sqrt(dt));
    for (int j = 1; j < m; j++)
    {
        Z[j] = (log(S[j] / S[j - 1]) - (r - sigma * sigma / 2.0) * dt) / (sigma * sqrt(dt));
    }
}

void Rescale(SamplePath &S, vector<double> &Z, double S0,
             double r, double sigma, double dt)
{
    int m = S.size();
    S[0] = S0 * exp((r - sigma * sigma / 2.0) * dt + sigma * sqrt(dt) * Z[0]);
    for (int j = 1; j < m; j++)
    {
        S[j] = S[j - 1] * exp((r - sigma * sigma / 2.0) * dt + sigma * sqrt(dt) * Z[j]);
    }
}

double PathDepOption::PriceByMC(BSModel Model, long N,
                                double epsilon)
{
    double H = 0.0, Hsq = 0.0, Hdelta = 0.0, Hrho = 0.0, Hvega = 0.0, Htheta = 0.0, Hgamma = 0.0;
    double S0 = Model.S0, r = Model.r, sigma = Model.sigma, dt = T / m;

    vector<double> Z(m);
    SamplePath S(m);

    for (long i = 0; i < N; i++)
    {
        Model.GenerateSamplePath(T, m, S);
        H = (i * H + Payoff(S)) / (i + 1.0);
        Hsq = (i * Hsq + pow(Payoff(S), 2.0)) / (i + 1.0);

        GetZ(Z, S, S0, r, sigma, dt);

        Rescale(S, Z, S0 * (1.0 + epsilon), r, sigma, dt);
        Hdelta = (i * Hdelta + Payoff(S)) / (i + 1.0);

        Rescale(S, Z, S0, r * (1.0 + epsilon), sigma, dt);
        Hrho = (i * Hrho + Payoff(S)) / (i + 1.0);

        Rescale(S, Z, S0, r, sigma * (1.0 + epsilon), dt);
        Hvega = (i * Hvega + Payoff(S)) / (i + 1.0);

        Rescale(S, Z, S0, r, sigma, dt * (1.0 + epsilon));
        Htheta = (i * Htheta + Payoff(S)) / (i + 1.0);

        Rescale(S, Z, S0 * (1.0 - epsilon), r, sigma, dt);
        Hgamma = (i * Hgamma + Payoff(S)) / (i + 1.0);
    }
    Price = exp(-r * T) * H;
    PricingError = exp(-r * T) * sqrt(Hsq - H * H) / sqrt(N - 1.0);

    delta = exp(-r * T) * (Hdelta - H) / (S0 * epsilon);
    vega = exp(-r * T) * (Hvega - H) / (sigma * epsilon);
    theta = -(exp(-r * T * (1.0 + epsilon)) * Htheta - exp(-r * T) * H) / (T * epsilon);

    gamma = exp(-r * T) * (Hdelta - 2.0 * H + Hgamma) / (S0 * S0 * epsilon * epsilon);

    rho = (exp(-r * (1.0 + epsilon) * T) * Hrho - exp(-r * T) * H) / (r * epsilon);

    return Price;
}

double ArthmAsianCall::Payoff(SamplePath &S)
{
    double Ave = 0.0;
    for (int k = 0; k < m; k++)
    {
        Ave = (k * Ave + S[k]) / (k + 1.0);
    }
    if (Ave < K)
        return 0.0;
    return Ave - K;
}

double EurCall::Payoff(SamplePath &S)
{
    if (S[0] < K)
        return 0.0;
    return S[0] - K;
}

double EurPut::Payoff(SamplePath &S)
{
    if (K < S[0])
        return 0.0;
    return K - S[0];
}