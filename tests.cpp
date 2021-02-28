#include "header.h"

double HiSq(std::vector<double> &est, std::vector<double> &th, long long len, double d) {
    double hisq = 0;
    for (int i = 0; i < d; ++i) {
        if (th[i] > 0) {
            hisq += (est[i] - th[i] * len) * (est[i] - th[i] * len) / (th[i] * len);
        }
    }
    return hisq;
}

specialDistGenTest::specialDistGenTest() {
    this->testE(2);
    this->testBern(0.3);
    this->testBin(10, 0.3);
    this->testGeom(0.3);
    this->testPois(2);
    this->testU(1, 10);
    this->testN(2, 1);
}

double GoldSteinChiSquare(double alpha, double n) {
    double d;
    if (alpha <= 0.999 && alpha >= 0.5) {
        d = 2.0637 * pow((log(1.0 / (1 - alpha)) - 0.16), 0.4274) - 1.5774;
    }
    else {
        d = -2.0637 * pow((log(1.0 / (alpha)) - 0.16), 0.4274) + 1.5774;
    }
    double Chi = 0;
    vector <double> A = { 1.0000886,0.4713941,0.0001348028,-0.008553069,0.00312558,-0.0008426812,0.00009780499 };
    vector <double> B = { -0.2237368,0.02607083,0.01128186,-0.01153761,0.005169654,0.00253001,-0.001450117 };
    vector <double> C = { -0.01513904,-0.008986007,0.02277679,-0.01323293,-0.006950356,0.001060438,0.001565326 };
    for (int i = 0; i <= 6; ++i) {
        Chi += pow(n, -(float)i / 2) * pow(d, i) * (A[i] + B[i] / n + C[i] / (n * n));
    }
    Chi = n * Chi * Chi * Chi;
    return Chi;
}

double specialDistGenTest::testBern(const double p) {

    std::vector<double> sampleS(size, 0);
    long d(floor(sqrt(sampleS.size()/5)));
    std::vector<double> sample(size, 0), probability(d, 0);
    for (long int i(0); i < size; ++i) {
        sampleS[i] = generatorS.generateBern(p);
    }

    std::vector <double> res(d+1);
    for (long int i(0); i < sample.size(); ++i) {
        res[sampleS[i]] += 1;
    }
    probability[0] = p;
    probability[1] = 1 - p;
    double V = HiSq(res, probability, sample.size(), d - 2);
    cout << "Статистика V равна: " << V << endl;
    cout << "Введите альфу, для которой хотите проверить гипотезу: ";
    double alpha = 0.95;
    double ChiQuant = GoldSteinChiSquare(alpha, d-2-1);
    cout << "Хи-квадрат для альфа равное: " << ChiQuant << endl;
    if (V < ChiQuant) {
        cout << "Гипотеза принимается" << endl;
    }
    else {
        cout << "Гипотеза отвергается" << endl;
    }
    return V;
}

double specialDistGenTest::testBin(const double n, const double p) {

    std::vector<double> sampleS(size, 0);
    long d(n);
    double max(0);

    std::vector<double> sample(size, 0), probability(d, 0);
    for (long int i(0); i < size; ++i) {
        sampleS[i] = generatorS.generateBin(n, p);
        if (sampleS[i] > max) {
            max = sampleS[i];
        }
    }

    std::vector <double> res(d+1);

    for (long int i(0); i < sampleS.size(); ++i) {
        res[sampleS[i]] += 1;
    }

    for (long int i(0); i < d-1; ++i) {
        probability[i] = pow(p, n - i)*pow(1 - p, i)*fact(n)/(fact(i)*fact(n - i));
    }
    double V = HiSq(res, probability, sample.size(), d - 2);
    cout << "Статистика V равна: " << V << endl;
    cout << "Введите альфу, для которой хотите проверить гипотезу: ";
    double alpha = 0.95;
    double ChiQuant = GoldSteinChiSquare(alpha, d-2-1);
    cout << "Хи-квадрат для альфа равное: " << ChiQuant << endl;
    if (V < ChiQuant) {
        cout << "Гипотеза принимается" << endl;
    }
    else {
        cout << "Гипотеза отвергается" << endl;
    }
    return V;
}


double specialDistGenTest::testGeom(const double p) {

        std::vector<double> sampleS(size, 0);
        long d(floor(sqrt(sampleS.size()/5)));
        double max(0);

        std::vector<double> sample(size, 0), probability(d, 0);
        for (long int i(0); i < size; ++i) {
            sampleS[i] = generatorS.generateGeom(p);
            if (sampleS[i] > max) {
                max = sampleS[i];
            }
        }

        std::vector <double> res(d+1);

        for (long int i(0); i < sampleS.size(); ++i) {
            res[sampleS[i]] += 1;
        }

        for (long int i(0); i < d-1; ++i) {
            probability[i] = pow(1 - p, i)*p;
        }
    double V = HiSq(res, probability, sample.size(), d - 2);
    cout << "Статистика V равна: " << V << endl;
    cout << "Введите альфу, для которой хотите проверить гипотезу: ";
    double alpha = 0.95;
    double ChiQuant = GoldSteinChiSquare(alpha, d-2-1);
    cout << "Хи-квадрат для альфа равное: " << ChiQuant << endl;
    if (V < ChiQuant) {
        cout << "Гипотеза принимается" << endl;
    }
    else {
        cout << "Гипотеза отвергается" << endl;
    }
    return V;
}

double specialDistGenTest::testPois(const double l) {

    std::vector<double> sampleS(size, 0);
    long d(floor(sqrt(sampleS.size()/5)));
    double max(0);

    std::vector<double> sample(size, 0), probability(d, 0);
    for (long int i(0); i < size; ++i) {
        sampleS[i] = generatorS.generatePois(l);
        if (sampleS[i] > max) {
            max = sampleS[i];
        }
    }

    std::vector <double> res(d+1);

    for (long int i(0); i < sampleS.size(); ++i) {
        res[sampleS[i]] += 1;
    }

    for (long int i(0); i < d-1; ++i) {
        probability[i] = pow(l, i)*exp(-l)/fact(i);
    }
    double V = HiSq(res, probability, sample.size(), 30);
    cout << "Статистика V равна: " << V << endl;
    cout << "Введите альфу, для которой хотите проверить гипотезу: ";
    double alpha = 0.95;
    double ChiQuant = GoldSteinChiSquare(alpha, 30-1);
    cout << "Хи-квадрат для альфа равное: " << ChiQuant << endl;
    if (V < ChiQuant) {
        cout << "Гипотеза принимается" << endl;
    }
    else {
        cout << "Гипотеза отвергается" << endl;
    }
    return V;
}


double probability(double arg, double mean_, double variance_) {
    double pi = 3.14159265358979323846264;
    return std::max(0.0, std::exp(-arg * arg / (2 * variance_)) / std::sqrt(variance_ * 2 * pi));
}

double distribution(double arg, double mean_, double variance_) {
    double integral = 0;
    int steps_number = 100;
    double step = (arg - mean_) / steps_number;
    for (int i = 0; i < steps_number; ++i) {
        integral += probability((i + 0.5) * step, mean_, variance_) * step;
    }
    return integral + 0.5;
}

double specialDistGenTest::testN(const double m, const double s) {
    std::vector<double> sampleS(size, 0);
    long d(floor(sqrt(sampleS.size()/5)));
    double max(0), min(0);

    std::vector<double> sample(size, 0), probability(d + 10, 0);
    for (long int i(0); i < size; ++i) {
        sampleS[i] = generatorS.generateN(m, s);
        if (sampleS[i] > max) {
            max = sampleS[i];
        }
        if (sampleS[i] < min) {
            min = sampleS[i];
        }
    }
    for (long int i(0); i < sample.size(); ++i) {
        sample[i] = floor((d*(sampleS[i]- min))/(max - min));
        sampleS[i] = floor((d*(sampleS[i]- min))/(max - min))*(max - min)/d;
    }
    std::vector <double> res(d+10, 0);

    for (long int i(0); i < sample.size(); ++i) {
        res[sample[i]] += 1;
    }
    sort(sample.begin(), sample.end());
    std::unique(sample.begin(), sample.end());
    int z(0);
    for (double i(min); i < max; i += (max - min)/d) {
        probability[z] = distribution(i + (max - min)/d, m, s) - distribution(i, m, s);
        ++z;
    }
    double V = HiSq(res, probability, sample.size(), 30);
    cout << "Статистика V равна: " << V << endl;
    cout << "Введите альфу, для которой хотите проверить гипотезу: ";
    double alpha = 0.95;
    double ChiQuant = GoldSteinChiSquare(alpha, 30-1);
    cout << "Хи-квадрат для альфа равное: " << ChiQuant << endl;
    if (V < ChiQuant) {
        cout << "Гипотеза принимается" << endl;
    }
    else {
        cout << "Гипотеза отвергается" << endl;
    }
    return V;
}

double specialDistGenTest::testU(const double a, const double b) {
    std::vector<double> sampleS(size, 0);
    long d(floor(sqrt(sampleS.size()/5)));

    std::vector<double> sample(size, 0), probability(d, 0);
    for (long int i(0); i < size; ++i) {
        sampleS[i] = generatorS.generateU(a, b);
    }
    for (long int i(0); i < sample.size(); ++i) {
        sampleS[i] = floor(d*((sampleS[i] - a)/(b - a)));
    }

    std::vector <double> res(d+1, 0);
    for (long int i(0); i < sample.size(); ++i) {
        res[sampleS[i]] += 1;
    }
    std::vector<double> p(d, (double)1 / (d));
    double V = HiSq(res, probability, sample.size(), d-2);
    cout << "Статистика V равна: " << V << endl;
    cout << "Введите альфу, для которой хотите проверить гипотезу: ";
    double alpha = 0.95;
    double ChiQuant = GoldSteinChiSquare(alpha, d-2-1);
    cout << "Хи-квадрат для альфа равное: " << ChiQuant << endl;
    if (V < ChiQuant) {
        cout << "Гипотеза принимается" << endl;
    }
    else {
        cout << "Гипотеза отвергается" << endl;
    }
    return V;
}

double specialDistGenTest::testE(double l) {
    std::vector<double> sampleS(size, 0);
    long d(floor(sqrt(sampleS.size()/5)));
    double max(0);

    std::vector<double> sample(size, 0), probability(d, 0);
    for (long int i(0); i < size; ++i) {
        sampleS[i] = generatorS.generateE(1);
        if (sampleS[i] > max) {
            max = sampleS[i];
        }
    }
    for (long int i(0); i < sample.size(); ++i) {
        sample[i] = floor(d*sampleS[i]/max);
        sampleS[i] = floor(d*sampleS[i]/max)*max/d;
    }

    std::vector <double> res(d+1, 0);

    for (long int i(0); i < sample.size(); ++i) {
        res[sample[i]] += 1;
    }
    sort(sampleS.begin(), sampleS.end());
    std::unique(sampleS.begin(), sampleS.end());
    for (long int i(0); i < d-1; ++i) {
        probability[i] = (1-exp(-1*sampleS[i+1])) - (1 - exp(-1*sampleS[i]));
    }
    double V = HiSq(res, probability, sample.size(), 30);
    cout << "Статистика V равна: " << V << endl;
    cout << "Введите альфу, для которой хотите проверить гипотезу: ";
    double alpha = 0.95;
    double ChiQuant = GoldSteinChiSquare(alpha, 30-1);
    cout << "Хи-квадрат для альфа равное: " << ChiQuant << endl;
    if (V < ChiQuant) {
        cout << "Гипотеза принимается" << endl;
    }
    else {
        cout << "Гипотеза отвергается" << endl;
    }
    return V;
}

