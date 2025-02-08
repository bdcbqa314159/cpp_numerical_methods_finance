#include <iostream>

void bubble_sort(int a[], int size);

int main()
{
    const int size = 7;
    int a[size] = {10, 4, 56, 2, 3, 6, 8};

    for (int i = 0; i < 7; ++i)
    {
        std::cout << a[i] << " ";
    }

    std::cout << "\n";
    std::cout << "sorting...\n";
    bubble_sort(a, size);

    for (int i = 0; i < 7; ++i)
    {
        std::cout << a[i] << " ";
    }

    std::cout << "\n";

    return 0;
}

void bubble_sort(int a[], int size)
{
    for (int i = 0; i < size - 1; ++i)
    {
        for (int j = 0; j < size - i; ++j)
        {
            if (a[j] > a[j + 1])
            {
                int temp = a[j];
                a[j] = a[j + 1];
                a[j + 1] = temp;
            }
        }
    }
}