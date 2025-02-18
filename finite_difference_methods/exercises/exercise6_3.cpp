#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace std;

const double pi = 4.0 * atan(1.0);

typedef vector<double> Vector;
typedef std::vector<double> SamplePath;

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

class ExplicitMethod : public FDMethod
{
public:
    ExplicitMethod(ParabPDE *PtrPDE_, int imax_, int jmax_)
        : FDMethod(PtrPDE_, imax_, jmax_) {}

    double A(int i, int j)
    {
        return dt * (b(i, j) / 2.0 - a(i, j) / dx) / dx;
    }
    double B(int i, int j)
    {
        return 1.0 - dt * c(i, j) + 2.0 * dt * a(i, j) / (dx * dx);
    }
    double C(int i, int j)
    {
        return -dt * (b(i, j) / 2.0 + a(i, j) / dx) / dx;
    }
    double D(int i, int j)
    {
        return -dt * d(i, j);
    }

    void SolvePDE();
};

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

double x(FDMethod *Method, double t, double S)
{
    return (Method->v(t, S + Method->dx) - Method->v(t, S - Method->dx)) / (2.0 * Method->dx);
}

double y(FDMethod *Method, double t, double S)
{
    return Method->v(t, S) - x(Method, t, S) * S;
}

int main()
{
    double S0 = 100.0, r = 0.05, sigma = 0.2;
    BSModel Model(S0, r, sigma);

    double K = 100.0, T = 1. / 12., zl = 0.0, zu = 2.0 * S0;
    Put EuropeanPut(K, T, zl, zu);

    BSEq BSPDE(&Model, &EuropeanPut);

    int imax = 3000, jmax = 1000;
    ExplicitMethod Method(&BSPDE, imax, jmax);

    Method.SolvePDE();

    double t = 0.0, S = 100.0;

    cout << "S = " << S << ", t = " << t << endl
         << endl;
    cout << "Price = " << Method.v(t, S) << endl
         << endl;
    cout << "Strategy:" << endl
         << "  number of stocks         = " << x(&Method, t, S) << endl
         << "  money invested risk free = " << y(&Method, t, S) << endl
         << endl;

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

void ExplicitMethod::SolvePDE()
{
    for (int j = 0; j <= jmax; j++)
        V[imax][j] = f(j);
    for (int i = imax; i > 0; i--)
    {
        V[i - 1][0] = fl(i - 1);
        V[i - 1][jmax] = fu(i - 1);
        for (int j = 1; j < jmax; j++)
        {
            V[i - 1][j] = A(i, j) * V[i][j - 1] + B(i, j) * V[i][j] + C(i, j) * V[i][j + 1] + D(i, j);
        }
    }
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
