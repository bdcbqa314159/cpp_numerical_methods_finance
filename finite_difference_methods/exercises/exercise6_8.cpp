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

class HeatEq : public ParabPDE
{
public:
    BSModel *PtrModel;
    Option *PtrOption;
    HeatEq(BSModel *PtrModel_, Option *PtrOption_);

    double a(double t, double x) { return -0.5; }
    double b(double t, double x) { return 0.0; }
    double c(double t, double x) { return 0.0; }
    double d(double t, double x) { return 0.0; }

    double f(double x);
    double fl(double t);
    double fu(double t);

    double Z(double t, double x);
    double V(double t, double u);
    double X(double t, double z);
    double U(double t, double v);
};

class LCP
{
public:
    ParabPDE *PtrPDE;
    virtual double g(double t, double x) = 0;
};

class HeatEqLCP : public LCP, public HeatEq
{
public:
    HeatEqLCP(BSModel *PtrModel, Option *PtrOption)
        : HeatEq(PtrModel, PtrOption) { PtrPDE = this; }

    double g(double t, double x)
    {
        return V(t, PtrOption->Payoff(Z(t, x)));
    }
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

class FDLCP
{
public:
    LCP *PtrLCP;
    FDMethod *PtrFDMethod;
    double g(int i, int j)
    {
        return PtrLCP->g(PtrFDMethod->t(i), PtrFDMethod->x(j));
    }
};

class ExplicitLCP : public ExplicitMethod, public FDLCP
{
public:
    ExplicitLCP(LCP *PtrLCP_, int imax_, int jmax_)
        : ExplicitMethod(PtrLCP_->PtrPDE, imax_, jmax_)
    {
        PtrLCP = PtrLCP_;
        PtrFDMethod = this;
    }

    void SolveLCP();
};

int main()
{
    double S0 = 100.0, r = 0.05, sigma = 0.2;
    BSModel Model(S0, r, sigma);

    double K = 100.0, T = 1. / 12., zl = 1.0, zu = 2.0 * S0;
    Put PutOption(K, T, zl, zu);

    HeatEqLCP HeatLCP(&Model, &PutOption);

    int imax = 3000, jmax = 1000;
    ExplicitLCP Method(&HeatLCP, imax, jmax);

    Method.SolveLCP();

    double t = 0.0;
    double z = S0;
    double x = HeatLCP.X(t, z);
    cout << "Am  Put Price = " << HeatLCP.U(t, Method.v(t, x)) << endl;

    Method.SolvePDE();
    cout << "Eur Put Price = " << HeatLCP.U(t, Method.v(t, x)) << endl;

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

HeatEq::HeatEq(BSModel *PtrModel_, Option *PtrOption_)
{
    PtrModel = PtrModel_;
    PtrOption = PtrOption_;
    T = PtrOption->T;
    xl = X(0.0, PtrOption->zl);
    xu = X(0.0, PtrOption->zu);
}

double HeatEq::f(double x)
{
    return V(T, PtrOption->Payoff(Z(T, x)));
}

double HeatEq::fl(double t)
{
    return V(t, PtrOption->LowerBdCond(PtrModel, t));
}

double HeatEq::fu(double t)
{
    return V(t, PtrOption->UpperBdCond(PtrModel, t));
}

double HeatEq::Z(double t, double x)
{
    double r = PtrModel->r;
    double sigma = PtrModel->sigma;
    double S0 = PtrModel->S0;
    return S0 * exp((r - 0.5 * sigma * sigma) * t + sigma * x);
}

double HeatEq::V(double t, double u)
{
    return exp(-PtrModel->r * t) * u;
}

double HeatEq::X(double t, double z)
{
    double r = PtrModel->r;
    double sigma = PtrModel->sigma;
    double S0 = PtrModel->S0;
    return (log(z / S0) - (r - 0.5 * sigma * sigma) * t) / sigma;
}

double HeatEq::U(double t, double v)
{
    return exp(PtrModel->r * t) * v;
}

void ExplicitLCP::SolveLCP()
{
    for (int j = 0; j <= jmax; j++)
    {
        V[imax][j] = f(j);
        if (V[imax][j] < g(imax, j))
            V[imax][j] = g(imax, j);
    }
    for (int i = imax; i > 0; i--)
    {
        V[i - 1][0] = fl(i - 1);
        V[i - 1][jmax] = fu(i - 1);
        for (int j = 1; j < jmax; j++)
        {
            V[i - 1][j] = A(i, j) * V[i][j - 1] + B(i, j) * V[i][j] + C(i, j) * V[i][j + 1] + D(i, j);
        }
        for (int j = 0; j <= jmax; j++)
        {
            if (V[i - 1][j] < g(i - 1, j))
                V[i - 1][j] = g(i - 1, j);
        }
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
