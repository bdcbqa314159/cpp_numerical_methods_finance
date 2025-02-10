#pragma once
#include "binomial_model.hpp"
#include "binomial_lattice.hpp"
#include <iostream>

class Option
{
private:
    int N{};

public:
    Option() = default;
    Option(int N) : N(N) {}
    virtual ~Option() = default;

    void set_N(int N_)
    {
        N = N_;
    }

    int get_N() const { return N; }

    virtual double payoff(double x) const = 0;
};

class EurOption : public virtual Option
{
public:
    double price_by_CRR(BinomialModel Model);
};

class AmOption : public virtual Option
{
public:
    double price_by_Snell(BinomialModel Model, BinomialLattice<double> &priceTree, BinomialLattice<bool> &stoppingTree);
};

class SingleStrikeOption : public EurOption
{

protected:
    double K = 100.;

public:
    SingleStrikeOption() = default;
    SingleStrikeOption(int N, double K) : Option(N), K(K) {}

    void set_K(double K_) { K = K_; }
};

class Call : public SingleStrikeOption
{
public:
    Call() = default;

    double payoff(double x) const override;
};

class Put : public SingleStrikeOption
{
public:
    Put() = default;

    double payoff(double x) const override;
};

class DigitalCall : SingleStrikeOption
{
public:
    DigitalCall() = default;

    double payoff(double x) const override;
};

class DigitalPut : SingleStrikeOption
{
public:
    DigitalPut() = default;

    double payoff(double x) const override;
};

class DoubleStrikeOption : public EurOption
{
protected:
    double K1 = 100., K2 = 101.;

public:
    DoubleStrikeOption() = default;
    DoubleStrikeOption(int N, double K1, double K2) : Option(N), K1(K1), K2(K2) {}

    void set_K(double K1_, double K2_)
    {
        K1 = K1_;
        K2 = K2_;
    }
};

class BullSpread : public DoubleStrikeOption
{

public:
    double payoff(double x) const override;
};

class BearSpread : public DoubleStrikeOption
{

public:
    double payoff(double z) const override;
};

class DoubleDigital : public DoubleStrikeOption
{

public:
    double payoff(double x) const override;
};

class Strangle : public DoubleStrikeOption
{

public:
    double payoff(double x) const override;
};

class Butterfly : public DoubleStrikeOption
{

public:
    double payoff(double x) const override;
};