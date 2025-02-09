#include <iostream>

class DefInt
{

private:
    double a = 0., b = 1.;

public:
    virtual double function(double x) = 0;
    DefInt() = default;
    virtual ~DefInt() = default;

    void setBoundaries(double a, double b)
    {
        this->a = a;
        this->b = b;
    }

    double ByTrapezoid(int N)
    {
        double h = (b - a) / N;
        double integral = 0.;

        for (int i = 1; i < N; ++i)
        {
            double x_i = a + i * h;
            integral += (2 * function(x_i));
        }
        integral += (function(a) + function(b));
        integral *= (0.5 * h);
        return integral;
    }

    double BySimpson(int N)
    {
        double h = (b - a) / N;
        double integral = 0;

        for (int n = 1; n < N; n++)
            integral += 4 * function(a + n * h - 0.5 * h) + 2 * function(a + n * h);

        integral += 4 * function(b - 0.5 * h) + function(b);
        integral *= h / 6;
        return integral;
    }
};

class MyIntegral : public DefInt
{
public:
    virtual double function(double x) override
    {
        return x * x;
    }
};

int main()
{
    double a = 0., b = 1.;
    MyIntegral myInt;
    myInt.setBoundaries(a, b);
    int N = 100;

    std::cout << "Integral of x*x on the interval [0,1]:\n";
    std::cout << myInt.ByTrapezoid(N) << "\n";
    std::cout << myInt.BySimpson(N) << "\n";

    return 0;
}

double f(double x)
{
    return x * x;
}