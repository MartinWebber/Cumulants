#pragma once
#include <complex>

class __declspec(dllexport) Cumulants {
private:
    std::complex <float> *cumulants;
public:
    Cumulants(std::complex <float> *cumulants);
    std::complex <float> *getCumulants();
    void setCumulants(std::complex <float> *cumulants);
    std::complex <float> epr(unsigned p, unsigned r, unsigned k, unsigned n, std::complex <float> *sPhase, std::complex <float> *sComp);
    std::complex <float> cpr(unsigned p, unsigned r, unsigned k, unsigned n, std::complex <float> *sPhase, std::complex <float> *sComp);
    void cprAll(unsigned k, unsigned n, std::complex <float> *sPhase, std::complex <float> *sComp);
};