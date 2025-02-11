#include <iostream>
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

class F1
{
public:
    double Value(double x) { return x * x - 2; }
    double Deriv(double x) { return 2 * x; }
} MyF1;

class F2
{
private:
    double a; // parameter
public:
    F2(double a_) { a = a_; }
    double Value(double x) { return x * x - a; }
    double Deriv(double x) { return 2 * x; }
} MyF2(3.0);

int main()
{
    double Acc = 0.001;
    double LEnd = 0.0, REnd = 2.0;
    double Tgt = 0.0;
    cout << "Root of F1 by bisect: "
         << SolveByBisect(MyF1, Tgt, LEnd, REnd, Acc)
         << endl;
    cout << "Root of F2 by bisect: "
         << SolveByBisect(MyF2, Tgt, LEnd, REnd, Acc)
         << endl;
    double Guess = 1.0;
    cout << "Root of F1 by Newton-Raphson: "
         << SolveByNR(MyF1, Tgt, Guess, Acc)
         << endl;
    cout << "Root of F2 by Newton-Raphson: "
         << SolveByNR(MyF2, Tgt, Guess, Acc)
         << endl;
    return 0;
}
