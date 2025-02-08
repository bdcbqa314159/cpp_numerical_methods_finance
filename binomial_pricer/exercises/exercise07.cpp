/*
set(exe1 exercise01)
add_executable(${exe1} ${exe1}.cpp)
*/

#include <iostream>

void bubble_sort(int a[], int size);
void interchange(int *a, int *b);

int main()
{
    int a = 1, b = 3;
    std::cout << "before:\n"
              << "a = " << a << "\nb = " << b << "\n";
    interchange(&a, &b);
    std::cout << "after:\n"
              << "a = " << a << "\nb = " << b << "\n";

    const int size = 7;
    int x[size] = {10, 4, 56, 2, 3, 6, 8};

    for (int i = 0; i < 7; ++i)
    {
        std::cout << x[i] << " ";
    }

    std::cout << "\n";
    std::cout << "sorting...\n";
    bubble_sort(x, size);

    for (int i = 0; i < 7; ++i)
    {
        std::cout << x[i] << " ";
    }

    std::cout << "\n";

    return 0;
}

void interchange(int *a, int *b)
{
    if (a && b)
    {
        int temp = *a;
        *a = *b;
        *b = temp;
    }
}

void bubble_sort(int a[], int size)
{
    for (int i = 0; i < size - 1; ++i)
    {
        for (int j = 0; j < size - i; ++j)
        {
            if (a[j] > a[j + 1])
            {
                interchange(&a[j], &a[j + 1]);
            }
        }
    }
}