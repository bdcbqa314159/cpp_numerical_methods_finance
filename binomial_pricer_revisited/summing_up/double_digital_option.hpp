#ifndef DoubDigitOpt_h
#define DoubDigitOpt_h

#include "options06.hpp"

class DoubDigitOpt : public EurOption
{
private:
    double K1; // parameter 1
    double K2; // parameter 2
public:
    int GetInputData();
    double Payoff(double z);
};

#endif
