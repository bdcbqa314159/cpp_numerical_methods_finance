/*
set(exe1 exercise01)
add_executable(${exe1} ${exe1}.cpp)
*/

#include <iostream>

void interchange(int &a, int &b)
{
    int temp = a;
    a = b;
    b = temp;
}

int main()
{
    int a = 1, b = 3;
    std::cout << "before:\n"
              << "a = " << a << "\nb = " << b << "\n";
    interchange(a, b);
    std::cout << "after:\n"
              << "a = " << a << "\nb = " << b << "\n";

    return 0;
}