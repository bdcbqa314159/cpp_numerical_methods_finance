#ifndef ExplicitLCP_h
#define ExplicitLCP_h

#include "lcp.hpp"
#include "fd_lcp.hpp"
#include "explicit_method.hpp"

class ExplicitLCP : public ExplicitMethod, public FDLCP
{
public:
    ExplicitLCP(LCP *PtrLCP_, int imax_, int jmax_)
        : ExplicitMethod(PtrLCP_->PtrPDE, imax_, jmax_)
    {
        PtrLCP = PtrLCP_;
        PtrFDMethod = this;
    }

    void SolveLCP();
};

#endif
