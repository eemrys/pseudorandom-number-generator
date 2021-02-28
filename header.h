#pragma once
#ifndef LWII_HEADER_H
#define LWII_HEADER_H
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <set>

using namespace std;

class uniformDistGen {
public:
    uniformDistGen() = default;
    ~uniformDistGen() = default;
    double random();
private:
      long long a{3141592621}, b{1}, m{static_cast<long long>(pow(10, 10))}, seed{0};
};

class specialDistGen {
public:
    specialDistGen();
    ~specialDistGen() = default;
    double generateBern(const double p);
    double generateBin(const double n, const double p);
    double generateGeom(const double p);
    double generateGeom(const double w, const double k, const double n);
    double generatePois(const double l);
    double generateU(const double a, const double b);
    double generateN(const double m, const double s);
    double generateE(const double l);
protected:
    uniformDistGen generator;
};

class specialDistGenTest : public specialDistGen {
public:
    specialDistGenTest();
    ~specialDistGenTest() = default;
    double testBern(const double p);
    double testBin(const double n, const double p);
    double testGeom(const double p);
    double testGeom(const double w, const double b, const double n);
    double testPois(const double l);
    double testU(const double a, const double b);
    double testN(const double m, const double v);
    double testE(double l);
protected:
    long long size{12000};
    specialDistGen generatorS;
};

long long inline fact(long long n) {
    if (n <= 1)
        return 1;
    else
        return n * fact(n - 1);
}

long long inline C(long long n, long long k) {
    return fact(n)/(fact(k)*fact(n - k));
}

#endif
