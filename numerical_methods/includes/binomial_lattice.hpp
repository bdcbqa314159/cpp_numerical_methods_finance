#pragma once

#include <iostream>
#include <iomanip>
#include <vector>

template <typename T>
class BinomialLattice
{
private:
    int N = 0;
    std::vector<std::vector<T>> lattice;

public:
    BinomialLattice() = default;

    void set_N(int N_)
    {
        N = N_;
        lattice.resize(N + 1);
        for (int n = 0; n <= N; ++n)
            lattice[n].resize(n + 1);
    }

    void set_node(int n, int i, T x)
    {
        lattice[n][i] = x;
    }

    T get_node(int n, int i)
    {
        return lattice[n][i];
    }

    void display() const
    {
        std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(3);
        for (int n = 0; n <= N; ++n)
        {
            for (int i = 0; i <= n; ++i)
                std::cout << std::setw(7) << get_node(n, i);
            std::cout << "\n";
        }
        std::cout << "\n";
    }
};