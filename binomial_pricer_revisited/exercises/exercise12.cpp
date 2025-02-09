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

    double BySimpson()
    {
        double h = (b - a) / 6;
        double mid_point = a + (b - a) * 0.5;

        double integral = h * (function(a) + 4 * function(mid_point) + function(b));
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
    std::cout << myIntegral.BySimpson() << "\n";

    return 0;
}

double f(double x)
{
    return x * x;
}