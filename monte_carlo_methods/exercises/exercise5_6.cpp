#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>

using namespace std;
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

class PathDepOption
{
public:
    double T, Price, PricingError, delta;
    int m;
    virtual double Payoff(SamplePath &S) = 0;
    double PriceByMC(BSModel Model, long N, double epsilon);
    double PriceByVarRedMC(BSModel Model, long N, PathDepOption &CVOption, double epsilon);
    virtual double PriceByBSFormula(BSModel Model) { return 0.0; }
    virtual double DeltaByBSFormula(BSModel Model) { return 0.0; }
};

class DifferenceOfOptions : public PathDepOption
{
public:
    PathDepOption *Ptr1;
    PathDepOption *Ptr2;
    DifferenceOfOptions(double T_, int m_,
                        PathDepOption *Ptr1_,
                        PathDepOption *Ptr2_)
    {
        T = T_;
        m = m_;
        Ptr1 = Ptr1_;
        Ptr2 = Ptr2_;
    }
    double Payoff(SamplePath &S)
    {
        return Ptr1->Payoff(S) - Ptr2->Payoff(S);
    }
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

class BarrierCall : public PathDepOption
{
public:
    double K, L;
    BarrierCall(double T_, double K_, double L_, int m_)
    {
        T = T_;
        K = K_;
        L = L_;
        m = m_;
    }
    double Payoff(SamplePath &S);
};

class EurCall : public PathDepOption
{
public:
    double T, K;
    EurCall(double T_, double K_)
    {
        T = T_;
        K = K_;
        m = 1;
    }
    EurCall(double T_, double K_, int m_)
    {
        T = T_;
        K = K_;
        m = m_;
    }
    double d_plus(double S0, double sigma, double r);
    double d_minus(double S0, double sigma, double r);
    double PriceByBSFormula(double S0,
                            double sigma, double r);
    double VegaByBSFormula(double S0,
                           double sigma, double r);
    double DeltaByBSFormula(double S0,
                            double sigma, double r);
    double PriceByBSFormula(BSModel Model);
    double DeltaByBSFormula(BSModel Model);
    double Payoff(SamplePath &S);
};

int main()
{
    double S0 = 100.0, r = 0.03, sigma = 0.2;
    BSModel Model(S0, r, sigma);

    double T = 1.0 / 12.0, K = 100.0, L = 120.0;
    int m = 30;

    BarrierCall Option(T, K, L, m);
    EurCall CVOption(T, K, m);

    long N = 30000;
    double epsilon = 0.001;
    Option.PriceByVarRedMC(Model, N, CVOption, epsilon);
    cout << "Barrier call price = " << Option.Price << endl
         << "delta = " << Option.delta << endl
         << "Error = " << Option.PricingError << endl
         << endl;

    Option.PriceByMC(Model, N, epsilon);
    cout << "Price by direct MC = " << Option.Price << endl
         << "delta    = " << Option.delta << endl
         << "MC Error = " << Option.PricingError << endl;

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
    for (int j = 0; j < m; j++)
        S[j] = x * S[j];
}

double PathDepOption::PriceByMC(BSModel Model, long N, double epsilon)
{
    double H = 0.0, Hsq = 0.0, Heps = 0.0;
    SamplePath S(m);
    for (long i = 0; i < N; i++)
    {
        Model.GenerateSamplePath(T, m, S);
        H = (i * H + Payoff(S)) / (i + 1.0);
        Hsq = (i * Hsq + pow(Payoff(S), 2.0)) / (i + 1.0);
        Rescale(S, 1.0 + epsilon);
        Heps = (i * Heps + Payoff(S)) / (i + 1.0);
    }
    Price = exp(-Model.r * T) * H;
    PricingError = exp(-Model.r * T) * sqrt(Hsq - H * H) / sqrt(N - 1.0);
    delta = exp(-Model.r * T) * (Heps - H) / (Model.S0 * epsilon);
    return Price;
}

double PathDepOption::PriceByVarRedMC(BSModel Model, long N, PathDepOption &CVOption, double epsilon)
{
    DifferenceOfOptions VarRedOpt(T, m, this, &CVOption);

    Price = VarRedOpt.PriceByMC(Model, N, epsilon) + CVOption.PriceByBSFormula(Model);

    delta = VarRedOpt.delta + CVOption.DeltaByBSFormula(Model);

    PricingError = VarRedOpt.PricingError;

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

double BarrierCall::Payoff(SamplePath &S)
{
    for (int k = 0; k < m; k++)
    {
        if (S[k] > L)
            return 0.0;
    }
    if (S[m - 1] < K)
        return 0.0;
    return S[m - 1] - K;
}

double N(double x)
{
    double gamma = 0.2316419;
    double a1 = 0.319381530;
    double a2 = -0.356563782;
    double a3 = 1.781477937;
    double a4 = -1.821255978;
    double a5 = 1.330274429;
    double pi = 4.0 * atan(1.0);
    double k = 1.0 / (1.0 + gamma * x);
    if (x >= 0.0)
    {
        return 1.0 - ((((a5 * k + a4) * k + a3) * k + a2) * k + a1) * k * exp(-x * x / 2.0) / sqrt(2.0 * pi);
    }
    else
        return 1.0 - N(-x);
}

double EurCall::d_plus(double S0, double sigma, double r)
{
    return (log(S0 / K) +
            (r + 0.5 * pow(sigma, 2.0)) * T) /
           (sigma * sqrt(T));
}

double EurCall::d_minus(double S0, double sigma, double r)
{
    return d_plus(S0, sigma, r) - sigma * sqrt(T);
}

double EurCall::PriceByBSFormula(double S0,
                                 double sigma, double r)
{
    return S0 * N(d_plus(S0, sigma, r)) - K * exp(-r * T) * N(d_minus(S0, sigma, r));
}

double EurCall::VegaByBSFormula(double S0,
                                double sigma, double r)
{
    double pi = 4.0 * atan(1.0);
    return S0 * exp(-d_plus(S0, sigma, r) * d_plus(S0, sigma, r) / 2) * sqrt(T) / sqrt(2.0 * pi);
}

double EurCall::DeltaByBSFormula(double S0, double sigma, double r)
{
    return N(d_plus(S0, sigma, r));
}

double EurCall::PriceByBSFormula(BSModel Model)
{
    return PriceByBSFormula(Model.S0, Model.sigma, Model.r);
}

double EurCall::DeltaByBSFormula(BSModel Model)
{
    return DeltaByBSFormula(Model.S0, Model.sigma, Model.r);
}

double EurCall::Payoff(SamplePath &S)
{
    if (S[m - 1] < K)
        return 0.0;
    return S[m - 1] - K;
}
