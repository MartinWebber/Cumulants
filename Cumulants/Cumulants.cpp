#include "pch.h"
#include <iostream>
#include "Cumulants.h"

Cumulants::Cumulants(std::complex <float>* cumulants) {
    this->cumulants = cumulants;
}
std::complex <float> *Cumulants::getCumulants() {
    return cumulants;
}
void Cumulants::setCumulants(std::complex <float>* cumulants) {
    this->cumulants = cumulants;
}
std::complex <float> Cumulants::epr(unsigned p, unsigned r, unsigned n, std::complex <float>* sPhase, std::complex <float>* sComp) {
    std::complex <float> result = 0;
    for (unsigned i = 0; i < n; i++)
        result += pow(sPhase[i], p) * pow(sComp[i], r);
    result /= n;
    return result;
}
std::complex <float> Cumulants::cpr(unsigned p, unsigned r, unsigned n, std::complex <float>* sPhase, std::complex <float>* sComp) {
    std::complex <float> result = 0;
    if (p == 2 && (r == 0 || r == 1)) {
        result = epr(p, r, n, sPhase, sComp);
    }
    else if (p == 4 && r == 0) {
        result = epr(4, 0, n, sPhase, sComp) - 3.0f * (std::complex <float>) pow(epr(2, 0, n, sPhase, sComp), 2);
    }
    else if (p == 4 && r == 1) {
        result = epr(4, 1, n, sPhase, sComp)
            - 4.0f * epr(3, 0, n, sPhase, sComp) * epr(1, 1, n, sPhase, sComp)
            - 6.0f * epr(2, 1, n, sPhase, sComp) * epr(2, 0, n, sPhase, sComp);
    }
    else if (p == 4 && r == 2) {
        result = epr(4, 2, n, sPhase, sComp)
            - epr(4, 0, n, sPhase, sComp) * epr(0, 2, n, sPhase, sComp)
            - 8.0f * epr(3, 1, n, sPhase, sComp) * epr(1, 1, n, sPhase, sComp)
            - 4.0f * epr(3, 0, n, sPhase, sComp) * epr(1, 2, n, sPhase, sComp)
            - 6.0f * epr(2, 2, n, sPhase, sComp) * epr(2, 0, n, sPhase, sComp)
            - 6.0f * (std::complex <float>) pow(epr(2, 1, n, sPhase, sComp), 2)
            + 6.0f * (std::complex <float>) pow(epr(2, 0, n, sPhase, sComp), 2) * epr(0, 2, n, sPhase, sComp)
            + 24.0f * epr(2, 0, n, sPhase, sComp) * (std::complex <float>) pow(epr(1, 1, n, sPhase, sComp), 2);
    }
    else if (p == 6 && r == 0) {
        result = epr(6, 0, n, sPhase, sComp)
            - 15.0f * epr(4, 0, n, sPhase, sComp) * epr(2, 0, n, sPhase, sComp)
            - 10.0f * (std::complex <float>) pow(epr(3, 0, n, sPhase, sComp), 2)
            + 30.0f * (std::complex <float>) pow(epr(2, 0, n, sPhase, sComp), 2) * epr(1, 1, n, sPhase, sComp);
    }
    else if (p == 6 && r == 1) {
        result = epr(6, 1, n, sPhase, sComp)
            - 15.0f * epr(4, 1, n, sPhase, sComp) * epr(2, 0, n, sPhase, sComp)
            - 6.0f * epr(5, 0, n, sPhase, sComp) * epr(1, 1, n, sPhase, sComp)
            - 20.0f * epr(3, 1, n, sPhase, sComp) * epr(3, 0, n, sPhase, sComp)
            - 15.0f * epr(4, 0, n, sPhase, sComp) * epr(2, 1, n, sPhase, sComp)
            + 90.0f * epr(2, 1, n, sPhase, sComp) * (std::complex <float>) pow(epr(2, 0, n, sPhase, sComp), 2)
            + 120.0f * epr(3, 0, n, sPhase, sComp) * epr(2, 0, n, sPhase, sComp) * epr(1, 1, n, sPhase, sComp);
    }
    else if (p == 7 && r == 1) {
        result = epr(7, 1, n, sPhase, sComp)
            - 21.0f * epr(5, 1, n, sPhase, sComp) * epr(2, 0, n, sPhase, sComp)
            - 7.0f * epr(1, 1, n, sPhase, sComp) * epr(6, 0, n, sPhase, sComp)
            - 35.0f * epr(4, 1, n, sPhase, sComp) * epr(3, 0, n, sPhase, sComp)
            - 21.0f * epr(2, 1, n, sPhase, sComp) * epr(5, 0, n, sPhase, sComp)
            - 35.0f * epr(4, 0, n, sPhase, sComp) * epr(3, 1, n, sPhase, sComp)
            + 210.0f * epr(4, 0, n, sPhase, sComp) * epr(2, 0, n, sPhase, sComp) * epr(1, 1, n, sPhase, sComp)
            + 210.0f * (std::complex <float>) pow(epr(2, 0, n, sPhase, sComp), 2) * epr(3, 1, n, sPhase, sComp)
            + 420.0f * epr(3, 0, n, sPhase, sComp) * epr(2, 1, n, sPhase, sComp) * epr(2, 0, n, sPhase, sComp)
            + 140.0f * (std::complex <float>) pow(epr(3, 0, n, sPhase, sComp), 2) * epr(1, 1, n, sPhase, sComp)
            - 630.0f * (std::complex <float>) pow(epr(2, 0, n, sPhase, sComp), 3) * epr(1, 1, n, sPhase, sComp);
    }
    else if (p == 8 && r == 0) {
        result = epr(8, 0, n, sPhase, sComp)
            - 28.0f * epr(6, 0, n, sPhase, sComp) * epr(2, 0, n, sPhase, sComp)
            - 56.0f * epr(5, 0, n, sPhase, sComp) * epr(3, 0, n, sPhase, sComp)
            - 35.0f * (std::complex <float>) pow(epr(4, 0, n, sPhase, sComp), 2)
            + 420.0f * epr(4, 0, n, sPhase, sComp) * (std::complex <float>) pow(epr(2, 0, n, sPhase, sComp), 2)
            + 560.0f * (std::complex <float>) pow(epr(3, 0, n, sPhase, sComp), 2) * epr(2, 0, n, sPhase, sComp)
            - 630.0f * (std::complex <float>) pow(epr(2, 0, n, sPhase, sComp), 4);
    }
    else std::cout << "Error";
    return result;
}
void Cumulants::cprAll(unsigned unsigned n, std::complex <float>* sPhase, std::complex <float>* sComp) {
    std::complex <float> e20 = epr(2, 0, n, sPhase, sComp),
        e21 = epr(2, 1, n, sPhase, sComp),
        e40 = epr(4, 0, n, sPhase, sComp),
        e41 = epr(4, 1, n, sPhase, sComp),
        e30 = epr(3, 0, n, sPhase, sComp),
        e11 = epr(1, 1, n, sPhase, sComp),
        e42 = epr(4, 2, n, sPhase, sComp),
        e02 = epr(0, 2, n, sPhase, sComp),
        e31 = epr(3, 1, n, sPhase, sComp),
        e12 = epr(1, 2, n, sPhase, sComp),
        e22 = epr(2, 2, n, sPhase, sComp),
        e60 = epr(6, 0, n, sPhase, sComp),
        e61 = epr(6, 1, n, sPhase, sComp),
        e50 = epr(5, 0, n, sPhase, sComp),
        e71 = epr(7, 1, n, sPhase, sComp),
        e51 = epr(5, 1, n, sPhase, sComp),
        e80 = epr(8, 0, n, sPhase, sComp);
    cumulants[0] = e20;
    cumulants[1] = e21;
    cumulants[2] = e40 - 3.0f * (std::complex <float>) pow(e20, 2);
    cumulants[3] = e41 - 4.0f * e30 * e11 - 6.0f * e21 * e20;
    cumulants[4] = e42 - e40 * e02 - 8.0f * e31 * e11 - 4.0f * e30 * e12 - 6.0f * e22 * e20 - 6.0f * (std::complex <float>) pow(e21, 2) + 6.0f * (std::complex <float>) pow(e20, 2) * e02 + 24.0f * e20 * (std::complex <float>) pow(e11, 2);
    cumulants[5] = e60 - 15.0f * e40 * e20 - 10.0f * (std::complex <float>) pow(e30, 2) + 30.0f * (std::complex <float>) pow(e20, 2) * e11;
    cumulants[6] = e61 - 15.0f * e41 * e20 - 6.0f * e50 * e11 - 20.0f * e31 * e30 - 15.0f * e40 * e21 + 90.0f * e21 * (std::complex <float>) pow(e20, 2) + 120.0f * e30 * e20 * e11;
    cumulants[7] = e71 - 21.0f * e51 * e20 - 7.0f * e11 * e60 - 35.0f * e41 * e30 - 21.0f * e21 * e50 - 35.0f * e40 * e31 + 210.0f * e40 * e20 * e11 + 210.0f * (std::complex <float>) pow(e20, 2) * e31 + 420.0f * e30 * e21 * e20 + 140.0f * (std::complex <float>) pow(e30, 2) * e11 - 630.0f * (std::complex <float>) pow(e20, 3) * e11;
    cumulants[8] = e80 - 28.0f * e60 * e20 - 56.0f * e50 * e30 - 35.0f * (std::complex <float>) pow(e40, 2) + 420.0f * e40 * (std::complex <float>) pow(e20, 2) + 560.0f * (std::complex <float>) pow(e30, 2) * e20 - 630.0f * (std::complex <float>) pow(e20, 4);
}