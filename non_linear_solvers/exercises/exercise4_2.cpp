#include <iostream>
#include <cassert>
using namespace std;

template <typename Function>
double SolveByBisect(Function &Fct,
                     double Tgt, double LEnd, double REnd, double Acc)
{
    double left = LEnd, right = REnd, mid = (left + right) / 2;
    double y_left = Fct.Value(left) - Tgt, y_mid = Fct.Value(mid) - Tgt;
    while (mid - left > Acc)
    {
        if ((y_left > 0 && y_mid > 0) || (y_left < 0 && y_mid < 0))
        {
            left = mid;
            y_left = y_mid;
        }
        else
            right = mid;
        mid = (left + right) / 2;
        y_mid = Fct.Value(mid) - Tgt;
    }
    return mid;
}

template <typename Function>
double SolveByNR(Function &Fct,
                 double Tgt, double Guess, double Acc)
{
    double x_prev = Guess;
    double x_next = x_prev - (Fct.Value(x_prev) - Tgt) / Fct.Deriv(x_prev);
    while (x_next - x_prev > Acc || x_prev - x_next > Acc)
    {
        x_prev = x_next;
        x_next = x_prev - (Fct.Value(x_prev) - Tgt) / Fct.Deriv(x_prev);
    }
    return x_next;
}

class BondYield
{
public:
    std::vector<double> coupons, times;
    double face_value{};

    BondYield() = default;
    BondYield(const std::vector<double> &coupons, const std::vector<double> &times, double face_value) : coupons(coupons), times(times), face_value(face_value)
    {
        assert(coupons.size() == times.size());
    }

    double Value(double x)
    {
        double out{};
        for (size_t i = 0; i < coupons.size(); ++i)
        {
            out += coupons[i] * std::exp(-x * times[i]);
        }
        out += face_value * std::exp(-times[times.size() - 1] * x);
        return out;
    }

    double Deriv(double x)
    {
        double out{};
        for (size_t i = 0; i < coupons.size(); ++i)
        {
            out += (-times[i] * coupons[i] * std::exp(-x * times[i]));
        }
        out += (-times[times.size() - 1] * face_value * std::exp(-times[times.size() - 1] * x));
        return out;
    }
};

int main()
{
    double F = 100.0; // face value
    double T = 3.0;   // maturity time
    vector<double> C; // coupons
    C.push_back(1.2);
    C.push_back(1.2);
    C.push_back(1.2);
    vector<double> t; // coupon times
    t.push_back(1.0);
    t.push_back(2.0);
    t.push_back(3.0);

    BondYield MyBond(C, t, F);

    double P = 98.56;
    double Acc = 0.0001;
    double y;

    cout << "F = " << F << endl;
    cout << "T = " << F << endl;
    cout << "coupons: " << endl;
    for (unsigned int n = 0; n < C.size(); n++)
        cout << "C" << n << " = " << C[n] << " " << endl;
    cout << "tenors: " << endl;
    for (unsigned int n = 0; n < t.size(); n++)
        cout << "T" << n << " = " << t[n] << " " << endl;
    cout << "P = " << P << endl
         << endl;

    double LEnd = 0.0;
    double REnd = 1.0;
    y = SolveByBisect(MyBond, P, LEnd, REnd, Acc);
    cout << "Yield by bisection method: " << y << endl;

    double Guess = 0.2;
    y = SolveByNR(MyBond, P, Guess, Acc);
    cout << "Yield by Newton-Raphson method: " << y << endl;

    return 0;
}
