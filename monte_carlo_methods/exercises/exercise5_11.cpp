#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>

#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace std;

typedef std::vector<double> Vector;
typedef std::vector<Vector> Matrix;
typedef std::vector<Vector> SamplePath;
const double pi = 4.0 * atan(1.0);

Vector operator*(const Matrix &C, const Vector &V);
Vector operator*(const double &a, const Vector &V);
Vector operator+(const double &a, const Vector &V);
Vector operator+(const Vector &V, const Vector &W);
Vector operator*(const Vector &V, const Vector &W);
Vector exp(const Vector &V);
double operator^(const Vector &V, const Vector &W);

class BSModel
{
public:
    Vector S0, sigma;
    Matrix C;
    double r;
    BSModel(Vector S0_, double r_, Matrix C_);
    void GenerateSamplePath(double T, int m, SamplePath &S);
};

class PathDepOption
{
public:
    double T, Price, PricingError;
    int m;
    virtual double Payoff(SamplePath &S) = 0;
    double PriceByMC(BSModel Model, long N);
    double PriceByVarRedMC(BSModel Model, long N, PathDepOption &CVOption);
    virtual double PriceByBSFormula(BSModel Model) { return 0.0; }
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

class EurBasketCall : public PathDepOption
{
public:
    double K;
    EurBasketCall(double T_, double K_)
    {
        T = T_;
        K = K_;
        m = 1;
    }
    double Payoff(SamplePath &S);
};

class SumOfCalls : public PathDepOption
{
public:
    Vector K;
    SumOfCalls(double T_, Vector K_)
    {
        T = T_;
        K = K_;
        m = 1;
    }
    double Payoff(SamplePath &S);
    double PriceByBSFormula(BSModel Model);
};

class EurCall
{
public:
    double T, K;
    EurCall(double T_, double K_)
    {
        T = T_;
        K = K_;
    }
    double d_plus(double S0, double sigma, double r);
    double d_minus(double S0, double sigma, double r);
    double PriceByBSFormula(double S0,
                            double sigma, double r);
    double VegaByBSFormula(double S0,
                           double sigma, double r);
};

int main()
{
    int d = 3;
    Vector S0(d);
    S0[0] = 40.0;
    S0[1] = 60.0;
    S0[2] = 100.0;
    double r = 0.03;

    Matrix C(d);
    for (int i = 0; i < d; i++)
        C[i].resize(d);

    C[0][0] = 0.1;
    C[0][1] = -0.1;
    C[0][2] = 0.0;
    C[1][0] = -0.1;
    C[1][1] = 0.2;
    C[1][2] = 0.0;
    C[2][0] = 0.0;
    C[2][1] = 0.0;
    C[2][2] = 0.3;
    BSModel Model(S0, r, C);

    double T = 1.0 / 12.0, K = 200.0;

    EurBasketCall Option(T, K);

    // initiating the control variate:
    double V = 0.0;
    for (int j = 0; j < d; j++)
        V = V + Model.S0[j];
    Vector Kd = (K / V) * Model.S0;
    SumOfCalls CVOption(T, Kd);

    long N = 30000;
    Option.PriceByVarRedMC(Model, N, CVOption);
    cout << "European Basket Call Price using Control Variates = "
         << Option.Price << endl
         << "Pricing Error = " << Option.PricingError << endl
         << endl;

    Option.PriceByMC(Model, N);
    cout << "European Basket Call Price using direct MC        = "
         << Option.Price << endl
         << "Pricing Error = " << Option.PricingError << endl
         << endl;

    return 0;
}

Vector operator*(const Matrix &C, const Vector &V)
{
    int d = C.size();
    Vector W(d);
    for (int j = 0; j < d; j++)
    {
        W[j] = 0.0;
        for (int l = 0; l < d; l++)
            W[j] = W[j] + C[j][l] * V[l];
    }
    return W;
}

Vector operator+(const Vector &V, const Vector &W)
{
    int d = V.size();
    Vector U(d);
    for (int j = 0; j < d; j++)
        U[j] = V[j] + W[j];
    return U;
}

Vector operator+(const double &a, const Vector &V)
{
    int d = V.size();
    Vector U(d);
    for (int j = 0; j < d; j++)
        U[j] = a + V[j];
    return U;
}

Vector operator*(const double &a, const Vector &V)
{
    int d = V.size();
    Vector U(d);
    for (int j = 0; j < d; j++)
        U[j] = a * V[j];
    return U;
}

Vector operator*(const Vector &V, const Vector &W)
{
    int d = V.size();
    Vector U(d);
    for (int j = 0; j < d; j++)
        U[j] = V[j] * W[j];
    return U;
}

Vector exp(const Vector &V)
{
    int d = V.size();
    Vector U(d);
    for (int j = 0; j < d; j++)
        U[j] = exp(V[j]);
    return U;
}

double operator^(const Vector &V, const Vector &W)
{
    double sum = 0.0;
    int d = V.size();
    for (int j = 0; j < d; j++)
        sum = sum + V[j] * W[j];
    return sum;
}

double Gauss()
{
    double U1 = (rand() + 1.0) / (RAND_MAX + 1.0);
    double U2 = (rand() + 1.0) / (RAND_MAX + 1.0);
    return sqrt(-2.0 * log(U1)) * cos(2.0 * pi * U2);
}

Vector Gauss(int d)
{
    Vector Z(d);
    for (int j = 0; j < d; j++)
        Z[j] = Gauss();
    return Z;
}

BSModel::BSModel(Vector S0_, double r_, Matrix C_)
{
    S0 = S0_;
    r = r_;
    C = C_;
    srand(time(NULL));
    int d = S0.size();
    sigma.resize(d);
    for (int j = 0; j < d; j++)
        sigma[j] = sqrt(C[j] ^ C[j]);
}

void BSModel::GenerateSamplePath(double T, int m, SamplePath &S)
{
    Vector St = S0;
    int d = S0.size();
    for (int k = 0; k < m; k++)
    {
        S[k] = St * exp((T / m) * (r + (-0.5) * sigma * sigma) + sqrt(T / m) * (C * Gauss(d)));
        St = S[k];
    }
}

double PathDepOption::PriceByMC(BSModel Model, long N)
{
    double H = 0.0, Hsq = 0.0;
    SamplePath S(m);
    for (long i = 0; i < N; i++)
    {
        Model.GenerateSamplePath(T, m, S);
        H = (i * H + Payoff(S)) / (i + 1.0);
        Hsq = (i * Hsq + pow(Payoff(S), 2.0)) / (i + 1.0);
    }
    Price = exp(-Model.r * T) * H;
    PricingError = exp(-Model.r * T) * sqrt(Hsq - H * H) / sqrt(N - 1.0);
    return Price;
}

double PathDepOption::PriceByVarRedMC(BSModel Model, long N, PathDepOption &CVOption)
{
    DifferenceOfOptions VarRedOpt(T, m, this, &CVOption);

    Price = VarRedOpt.PriceByMC(Model, N) + CVOption.PriceByBSFormula(Model);

    PricingError = VarRedOpt.PricingError;

    return Price;
}

double ArthmAsianCall::Payoff(SamplePath &S)
{
    double Ave = 0.0;
    int d = S[0].size();
    Vector one(d);
    for (int i = 0; i < d; i++)
        one[i] = 1.0;
    for (int k = 0; k < m; k++)
    {
        Ave = (k * Ave + (one ^ S[k])) / (k + 1.0);
    }
    if (Ave < K)
        return 0.0;
    return Ave - K;
}

double EurBasketCall::Payoff(SamplePath &S)
{
    double Sum = 0.0;
    int d = S[0].size();
    for (int i = 0; i < d; i++)
        Sum = Sum + S[0][i];
    if (Sum < K)
        return 0.0;
    return Sum - K;
}

double SumOfCalls::Payoff(SamplePath &S)
{
    int d = S[0].size();
    double Sum = 0.0;
    for (int j = 0; j < d; j++)
    {
        if (S[0][j] > K[j])
            Sum = Sum + S[0][j] - K[j];
    }
    return Sum;
}

double SumOfCalls::PriceByBSFormula(BSModel Model)
{
    int d = Model.S0.size();
    double Sum = 0.0;
    for (int j = 0; j < d; j++)
    {
        EurCall Option(T, K[j]);
        Sum = Sum + Option.PriceByBSFormula(Model.S0[j], Model.sigma[j], Model.r);
    }
    return Sum;
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
