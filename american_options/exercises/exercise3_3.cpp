#include <iostream>

template <typename T>
void interchange(T &a, T &b)
{
    T temp = a;
    a = b;
    b = temp;
}

int main()
{
    int a = 1, b = 2;
    char c = 'c', d = 'd';
    double e = 3.456, f = 678.12;

    std::cout << "Before:\n";
    std::cout << "a = " << a << " b = " << b << "\n";
    std::cout << "c = " << c << " d = " << d << "\n";
    std::cout << "e = " << e << " f = " << f << "\n";

    interchange<int>(a, b);
    interchange<char>(c, d);
    interchange<double>(e, f);

    std::cout << "After interchange with template:\n";
    std::cout << "a = " << a << " b = " << b << "\n";
    std::cout << "c = " << c << " d = " << d << "\n";
    std::cout << "e = " << e << " f = " << f << "\n";

    return 0;
}