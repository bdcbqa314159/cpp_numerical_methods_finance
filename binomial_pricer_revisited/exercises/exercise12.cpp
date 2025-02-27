#include <iostream>

class DefInt
{

private:
    double a = 0., b = 1.;
    double (*function)(double x) = nullptr;

public:
    DefInt() = default;
    DefInt(double a, double b, double (*f)(double x)) : a(a), b(b), function(f) {}
    void setFunction(double (*f)(double x))
    {
        function = f;
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

double f(double x);

int main()
{
    double a = 0., b = 1.;
    DefInt myIntegral(a, b, f);
    int N = 100;

    std::cout << "Integral of x*x on the interval [0,1]:\n";
    std::cout << myIntegral.ByTrapezoid(N) << "\n";
    std::cout << myIntegral.BySimpson(N) << "\n";

    return 0;
}

double f(double x)
{
    return x * x;
}