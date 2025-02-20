#ifndef FDLCP_h
#define FDLCP_h

#include "fd_method.hpp"
#include "lcp.hpp"

class FDLCP
{
public:
    LCP *PtrLCP;
    FDMethod *PtrFDMethod;
    double g(int i, int j)
    {
        return PtrLCP->g(PtrFDMethod->t(i), PtrFDMethod->x(j));
    }
};

#endif
