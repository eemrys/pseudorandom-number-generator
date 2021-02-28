#include "header.h"

double uniformDistGen::random() {
    seed = (seed * a + b) % m;
    double tmp(m);
    tmp = fabs((double)seed)/tmp;
    return tmp;
}

double specialDistGen::generateBern(const double p) {
    double cake = generator.random();
    return (cake > p) ? 1 : 0;
}

double specialDistGen::generateBin(const double n, const double p) {
    while (true) {
        int current = floor(generator.random() * n);
        double supportVal = generator.random();
        double c(0.9);
        if (c * supportVal <
            pow(p, n - current) * pow(1 - p, current) * fact(n) / (fact(current) * fact(n - current))) {
            return current;
        }
    }
}

double specialDistGen::generateGeom(const double p) {
    int counter(0);
    while (true) {
        double supportVal = this->generateBern(1 - p);

        if (supportVal == 1) {
            return counter;
        }
        ++counter;
    }
}

double specialDistGen::generateGeom(const double n, const double K, const double N) {
    double k = floor(generator.random()*min(n,K));
    double result = (1 - C(n, k+1)*C(N - n, K - k - 1)/C(N, K));
    return result;
}

double specialDistGen::generatePois(const double l) {
    double supportValue (exp(-l)), p(1), k(0);
    do {
        ++k;
        p = p*(generator.random());
    } while ( p > supportValue);
    return k - 1;
}

double specialDistGen::generateU(const double a, const double b) {
    return generator.random()*(b-a) + a;
}

double specialDistGen::generateN(const double m, const double s) {
    double first, second, distance;
    do {
        first = 2 * generator.random() - 1;
        second = 2 * generator.random() - 1;
        distance = first * first + second * second;
    } while (distance > 1.0 || distance == 0);
    return std::sqrt(s) * first * std::sqrt(-2 * std::log(distance) / distance) + m;
}

double specialDistGen::generateE(const double l) {
    double cake = generator.random();
    return (-(1/l)*log(1 - cake));
}

specialDistGen::specialDistGen() {
    uniformDistGen gen;
    this->generator = gen;
}
