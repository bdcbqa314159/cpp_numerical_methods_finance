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
    double T, Price;
    int m;
    Vector delta;
    double PriceByMC(BSModel Model, long N, double epsilon);
    virtual ~PathDepOption() {}
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
    int m = 30;
    ArthmAsianCall Option(T, K, m);

    long N = 30000;
    double epsilon = 0.001;
    cout << "Arithmetic Basket Call Price = "
         << Option.PriceByMC(Model, N, epsilon) << endl;
    for (int j = 0; j < d; j++)
    {
        cout << "delta[" << j << "] = " << Option.delta[j] << endl;
    }

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

void Rescale(SamplePath &S, double x, int j)
{
    int m = S.size();
    for (int k = 0; k < m; k++)
        S[k][j] = x * S[k][j];
}

double PathDepOption::PriceByMC(BSModel Model, long N, double epsilon)
{
    double H = 0.0;
    SamplePath S(m);

    int d = Model.S0.size();
    delta.resize(d);

    Vector Heps(d);
    for (int i = 0; i < d; i++)
        Heps[i] = 0.0;

    for (long i = 0; i < N; i++)
    {
        Model.GenerateSamplePath(T, m, S);
        H = (i * H + Payoff(S)) / (i + 1.0);

        for (int j = 0; j < d; j++)
        {
            Rescale(S, 1.0 + epsilon, j);
            Heps[j] = (i * Heps[j] + Payoff(S)) / (i + 1.0);
            if (j < d - 1)
                Rescale(S, 1.0 / (1.0 + epsilon), j);
        }
    }
    Price = exp(-Model.r * T) * H;
    for (int j = 0; j < d; j++)
        delta[j] = exp(-Model.r * T) * (Heps[j] - H) / (epsilon * Model.S0[j]);
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
