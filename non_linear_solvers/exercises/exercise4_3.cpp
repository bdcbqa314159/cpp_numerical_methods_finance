#include <iostream>

template <typename Fct>
class DefInt
{

private:
    double a = 0., b = 1.;
    Fct function;

public:
    DefInt() = default;
    DefInt(double a, double b, Fct function) : a(a), b(b), function(function) {}
    void setFunction(Fct f)
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
            integral += (2 * function.value(x_i));
        }
        integral += (function.value(a) + function.value(b));
        integral *= (0.5 * h);
        return integral;
    }

    double BySimpson(int N)
    {
        double h = (b - a) / N;
        double integral = 0;

        for (int n = 1; n < N; n++)
            integral += 4 * function.value(a + n * h - 0.5 * h) + 2 * function.value(a + n * h);

        integral += 4 * function.value(b - 0.5 * h) + function.value(b);
        integral *= h / 6;
        return integral;
    }
};

class Function
{
public:
    double value(double x) { return x * x; }
};

int main()
{
    double a = 0., b = 1.;
    Function my_fct;
    DefInt myIntegral(a, b, my_fct);
    int N = 100;

    std::cout << "Integral of x*x on the interval [0,1]:\n";
    std::cout << myIntegral.ByTrapezoid(N) << "\n";
    std::cout << myIntegral.BySimpson(N) << "\n";

    return 0;
}
