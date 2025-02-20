#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>

using namespace std;
typedef vector<double> Vector;
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

class Call : public Option
{
public:
    double K;
    Call(double K_, double T_, double zl_, double zu_)
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

class ParabPDE
{
public:
    double T, xl, xu;

    virtual double a(double t, double x) = 0;
    virtual double b(double t, double x) = 0;
    virtual double c(double t, double x) = 0;
    virtual double d(double t, double x) = 0;

    virtual double f(double x) = 0;
    virtual double fu(double t) = 0;
    virtual double fl(double t) = 0;
};

class BSEq : public ParabPDE
{
public:
    BSModel *PtrModel;
    Option *PtrOption;
    BSEq(BSModel *PtrModel_, Option *PtrOption_);

    double a(double t, double z);
    double b(double t, double z);
    double c(double t, double z);
    double d(double t, double z);

    double f(double z);
    double fl(double t);
    double fu(double t);
};

class FDMethod
{
public:
    ParabPDE *PtrPDE;
    int imax, jmax;
    double dx, dt;

    vector<Vector> V;

    FDMethod(ParabPDE *PtrPDE_, int imax_, int jmax_);

    double t(double i) { return dt * i; }
    double x(int j) { return PtrPDE->xl + dx * j; }

    double a(double i, int j) { return PtrPDE->a(t(i), x(j)); }
    double b(double i, int j) { return PtrPDE->b(t(i), x(j)); }
    double c(double i, int j) { return PtrPDE->c(t(i), x(j)); }
    double d(double i, int j) { return PtrPDE->d(t(i), x(j)); }

    double f(int j) { return PtrPDE->f(x(j)); }
    double fu(int i) { return PtrPDE->fu(t(i)); }
    double fl(int i) { return PtrPDE->fl(t(i)); }

    double v(double t, double x);
};

typedef std::vector<double> Vector;
typedef std::vector<Vector> Matrix;

Vector operator*(const Matrix &C, const Vector &V);
Vector operator*(const double &a, const Vector &V);
Vector operator+(const double &a, const Vector &V);
Vector operator+(const Vector &V, const Vector &W);
Vector operator*(const Vector &V, const Vector &W);
Vector exp(const Vector &V);
double operator^(const Vector &V, const Vector &W);

class ImplicitScheme : public FDMethod
{
public:
    ImplicitScheme(ParabPDE *PtrPDE_, int imax_, int jmax_)
        : FDMethod(PtrPDE_, imax_, jmax_) {}

    virtual double A(int i, int j) = 0;
    virtual double B(int i, int j) = 0;
    virtual double C(int i, int j) = 0;
    virtual double D(int i, int j) = 0;
    virtual double E(int i, int j) = 0;
    virtual double F(int i, int j) = 0;
    virtual double G(int i, int j) = 0;

    Vector w(int i);
    Vector A(int i, Vector q);

    Vector LUDecomposition(int i, Vector q);

    void SolvePDE();
};

class ImplicitMethod : public ImplicitScheme
{
public:
    ImplicitMethod(ParabPDE *PtrPDE_, int imax_, int jmax_)
        : ImplicitScheme(PtrPDE_, imax_, jmax_)
    {
    }

    double A(int i, int j) { return 0.0; }
    double B(int i, int j) { return 1.0; }
    double C(int i, int j) { return 0.0; }
    double D(int i, int j) { return -dt * d(i - 1, j); }
    double E(int i, int j) { return -dt * (b(i - 1, j) / 2.0 - a(i - 1, j) / dx) / dx; }
    double F(int i, int j) { return 1.0 + dt * c(i - 1, j) - 2.0 * dt * a(i - 1, j) / (dx * dx); }
    double G(int i, int j) { return dt * (b(i - 1, j) / 2.0 + a(i - 1, j) / dx) / dx; }
};

class WeightImplicit : public ImplicitScheme
{
public:
    double lambda;
    WeightImplicit(ParabPDE *PtrPDE_, int imax_, int jmax_, double lambda_)
        : ImplicitScheme(PtrPDE_, imax_, jmax_) { lambda = lambda_; }

    double A(int i, int j)
    {
        return (1.0 - lambda) * dt * (b(i - lambda, j) / 2.0 - a(i - lambda, j) / dx) / dx;
    }
    double B(int i, int j)
    {
        return 1.0 + (1.0 - lambda) * dt * (2.0 * a(i - lambda, j) / (dx * dx) - c(i - lambda, j));
    }
    double C(int i, int j)
    {
        return -(1.0 - lambda) * dt * (b(i - lambda, j) / 2.0 + a(i - lambda, j) / dx) / dx;
    }

    double D(int i, int j) { return -dt * d(i - lambda, j); }

    double E(int i, int j)
    {
        return -lambda * dt * (b(i - lambda, j) / 2.0 - a(i - lambda, j) / dx) / dx;
    }
    double F(int i, int j)
    {
        return 1.0 - lambda * dt * (2.0 * a(i - lambda, j) / (dx * dx) - c(i - lambda, j));
    }
    double G(int i, int j)
    {
        return lambda * dt * (b(i - lambda, j) / 2.0 + a(i - lambda, j) / dx) / dx;
    }
};

int main()
{
    double S0 = 100.0, r = 0.05, sigma = 0.2;
    BSModel Model(S0, r, sigma);

    double T = 1. / 12., K = 100.0, zl = 0.0, zu = 2.0 * S0;
    Put EuropeanPut(K, T, zl, zu);

    int imax = 200, jmax = 2000;

    BSEq BSPDE(&Model, &EuropeanPut);

    double lambda = 0.6;
    WeightImplicit Method(&BSPDE, imax, jmax, lambda);
    Method.SolvePDE();
    cout << "Price = " << Method.v(0.0, S0) << endl;

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

double Call::Payoff(double z)
{
    if (z < K)
        return 0.0;
    return z - K;
}

double Call::UpperBdCond(BSModel *PtrModel, double t)
{
    return zu - K * exp(-PtrModel->r * (T - t));
}

double Call::LowerBdCond(BSModel *PtrModel, double t)
{
    return zl;
}

BSEq::BSEq(BSModel *PtrModel_, Option *PtrOption_)
{
    PtrModel = PtrModel_;
    PtrOption = PtrOption_;
    T = PtrOption->T;
    xl = PtrOption->zl;
    xu = PtrOption->zu;
}

double BSEq::a(double t, double z)
{
    return -0.5 * pow(PtrModel->sigma * z, 2.0);
}

double BSEq::b(double t, double z)
{
    return -PtrModel->r * z;
}

double BSEq::c(double t, double z)
{
    return PtrModel->r;
}

double BSEq::d(double t, double z)
{
    return 0.0;
}

double BSEq::f(double z)
{
    return PtrOption->Payoff(z);
}

double BSEq::fl(double t)
{
    return PtrOption->LowerBdCond(PtrModel, t);
}

double BSEq::fu(double t)
{
    return PtrOption->UpperBdCond(PtrModel, t);
}

Vector ImplicitScheme::w(int i)
{
    Vector w(jmax + 1);
    w[1] = D(i, 1) + A(i, 1) * fl(i) - E(i, 1) * fl(i - 1);
    for (int j = 2; j < jmax - 1; j++)
        w[j] = D(i, j);
    w[jmax - 1] = D(i, jmax - 1) + C(i, jmax - 1) * fu(i) - G(i, jmax - 1) * fu(i - 1);
    return w;
}

Vector ImplicitScheme::A(int i, Vector q)
{
    Vector p(jmax + 1);
    p[1] = B(i, 1) * q[1] + C(i, 1) * q[2];
    for (int j = 2; j < jmax - 1; j++)
    {
        p[j] = A(i, j) * q[j - 1] + B(i, j) * q[j] + C(i, j) * q[j + 1];
    }
    p[jmax - 1] = A(i, jmax - 1) * q[jmax - 2] + B(i, jmax - 1) * q[jmax - 1];
    return p;
}

Vector ImplicitScheme::LUDecomposition(int i, Vector q)
{
    Vector p(jmax + 1), r(jmax + 1), y(jmax + 1);
    r[1] = F(i, 1);
    y[1] = q[1];
    for (int j = 2; j < jmax; j++)
    {
        r[j] = F(i, j) - E(i, j) * G(i, j - 1) / r[j - 1];
        y[j] = q[j] - E(i, j) * y[j - 1] / r[j - 1];
    }
    p[jmax - 1] = y[jmax - 1] / r[jmax - 1];
    for (int j = jmax - 2; j > 0; j--)
    {
        p[j] = (y[j] - G(i, j) * p[j + 1]) / r[j];
    }
    return p;
}

void ImplicitScheme::SolvePDE()
{
    for (int j = 0; j <= jmax; j++)
        V[imax][j] = f(j);
    for (int i = imax; i > 0; i--)
    {
        V[i - 1] = LUDecomposition(i, A(i, V[i]) + w(i));
        V[i - 1][0] = fl(i - 1);
        V[i - 1][jmax] = fu(i - 1);
    }
}

FDMethod::FDMethod(ParabPDE *PtrPDE_, int imax_, int jmax_)
{
    PtrPDE = PtrPDE_;
    imax = imax_;
    jmax = jmax_;
    dx = (PtrPDE->xu - PtrPDE->xl) / jmax;
    dt = PtrPDE->T / imax;
    V.resize(imax + 1);
    for (int i = 0; i <= imax; i++)
        V[i].resize(jmax + 1);
}

double FDMethod::v(double t, double x)
{
    int i = (int)(t / dt);
    int j = (int)((x - PtrPDE->xl) / dx);
    double l1 = (t - FDMethod::t(i)) / dt, l0 = 1.0 - l1;
    double w1 = (x - FDMethod::x(j)) / dx, w0 = 1.0 - w1;
    return l1 * w1 * V[i + 1][j + 1] + l1 * w0 * V[i + 1][j] + l0 * w1 * V[i][j + 1] + l0 * w0 * V[i][j];
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
