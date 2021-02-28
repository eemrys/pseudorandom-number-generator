#include <iostream>
#include <vector>

using namespace std;

class NumberGenerator {
private:
    int a, b;
    long long m;
    vector<long long> Sequence;
public:
    long long X0;
    NumberGenerator(long long m, int a, int b, long long X0) {
        this->m = m;
        this->a = a;
        this->b = b;
        this->X0 = X0;
        vector<long long> Sequence;
    }

    long long GenerateNextNumber(long long PreviousNumber){
        return ((a * PreviousNumber + b) % m);
    }

    double Divisor(double Number)
    {
        return (Number / m);
    }
};


class Criteria {
public:
    Criteria Checker();

    double FrequencyTest(NumberGenerator Generator, int d, long long length) {
        vector <long long> rCounter(d,0);
        rCounter[Generator.X0]+=1;
        long long Un = Generator.X0;
        int Yn;
        for (long long i = 1; i < length; ++i)
        {
            Un = Generator.GenerateNextNumber(Un);
            Yn = (int)floor(d*Generator.Divisor(Un));
            rCounter[Yn]++;
        }
        vector <double> Prob(d, 1/(d*1.0));
        double Sum = Statistic(rCounter, length, d, Prob);
        return Sum;
    }

    double SerialTest(NumberGenerator Generator, int d, long long length) {
        vector <vector <long long> > rCounter(d, vector<long long>(d, 0));
        long long Un0 = Generator.X0;
        long long Un1 = Generator.GenerateNextNumber(Un0);
        int Yn0 = (int)floor(d * Generator.Divisor(Un0));
        int Yn1 = (int)floor(d * Generator.Divisor(Un1));
        rCounter[Yn0][Yn1]++;
        for (long long i = 1; i < length; ++i) {
            Un0 = Generator.GenerateNextNumber(Un1);
            Un1 = Generator.GenerateNextNumber(Un0);
            Yn0 = (int)floor(d * Generator.Divisor(Un0));
            Yn1 = (int)floor(d * Generator.Divisor(Un1));
            rCounter[Yn0][Yn1]++;
        }
        vector <long long> Counter(d*d, 0);
        double Sum;
        for (int i = 0; i < d; ++i)
        {
            for (int j = 0; j < d; ++j)
            {
                Counter[i * d + j] = rCounter[i][j];
            }
        }
        vector <double> Prob(d * d, 1 / (d * d * 1.0));
        Sum = Statistic(Counter, length, d * d, Prob);
        return Sum;
    }

    double CouponTest(NumberGenerator Generator, int d, long long length) {
        vector <long long> rCounter(6, 0);
        long long Un0 = Generator.X0;
        int Yn0 = (int)floor(d * Generator.Divisor(Un0));
        long long TempValueSequence = Un0;
        int TempValue = Yn0;
        long long SeqLength = 1;
        long long Un;
        int Yn;
        for (long long i = 1; i < length; ++i) {
            Un = Generator.GenerateNextNumber(TempValueSequence);
            Yn = (int)floor(d * Generator.Divisor(Un));
            if (Yn == TempValue) {
                SeqLength += 1;
            }
            else {
                if (SeqLength > 5) {
                    rCounter[5]++;
                }
                else {
                    rCounter[SeqLength - 1]++;
                }
                SeqLength = 1;
            }
            TempValue = Yn;
            TempValueSequence = Un;
        }
        if (SeqLength > 5) {
            rCounter[5]++;
        }
        else {
            rCounter[SeqLength - 1]++;
        }
        vector<double> Prob(6, 0);
        double TempMulti = (d-1);
        double SumProb = 0;
        long long AppearSum = 0; // cумма количества найденных последовательностей
        for (int i = 0; i < 5; ++i) {
            TempMulti /= (d * 1.0);
            Prob[i] = TempMulti;
            SumProb += Prob[i];
            AppearSum += rCounter[i];
        }
        Prob[5] = 1.0 - SumProb;
        double Sum = Statistic(rCounter, AppearSum, 6, Prob);
        return Sum;
    }

    double PokerTest(NumberGenerator Generator, int d, long long length) {
        vector<long long> rCounter(5, 0);
        long long Un0 = Generator.X0;
        long long Un1 = Generator.GenerateNextNumber(Un0);
        long long Un2 = Generator.GenerateNextNumber(Un1);
        long long Un3 = Generator.GenerateNextNumber(Un2);
        long long Un4 = Generator.GenerateNextNumber(Un3);
        long long NextValue = Un4;
        int Yn0 = (int)floor(d * Generator.Divisor(Un0));
        int Yn1 = (int)floor(d * Generator.Divisor(Un1));
        int Yn2 = (int)floor(d * Generator.Divisor(Un2));
        int Yn3 = (int)floor(d * Generator.Divisor(Un3));
        int Yn4 = (int)floor(d * Generator.Divisor(Un4));
        rCounter[Comparator(Yn0, Yn1, Yn2, Yn3, Yn4)]++;
        for (int i = 1; i < length / 5; ++i) {
            Un0 = Generator.GenerateNextNumber(NextValue);
            Un1 = Generator.GenerateNextNumber(Un0);
            Un2 = Generator.GenerateNextNumber(Un1);
            Un3 = Generator.GenerateNextNumber(Un2);
            Un4 = Generator.GenerateNextNumber(Un3);
            NextValue = Un4;
            Yn0 = (int)floor(d * Generator.Divisor(Un0));
            Yn1 = (int)floor(d * Generator.Divisor(Un1));
            Yn2 = (int)floor(d * Generator.Divisor(Un2));
            Yn3 = (int)floor(d * Generator.Divisor(Un3));
            Yn4 = (int)floor(d * Generator.Divisor(Un4));
            rCounter[Comparator(Yn0, Yn1, Yn2, Yn3, Yn4)]++;
        }
        vector<double> Prob { d*1 / (d*d*d*d*d * 1.0), d*(d-1)*15 / (d*d*d*d*d * 1.0), d*(d-1)*(d-2)*25 / (d*d*d*d*d * 1.0), d*(d-1)*(d-2)*(d-3)*10 / (d*d*d*d*d * 1.0), d*(d-1)*(d-2)*(d-3)*(d-4)*1 / (d*d*d*d*d * 1.0)};
        double Sum = Statistic(rCounter, length / 5, 5, Prob);
        return Sum;
    }

    double Statistic(vector <long long> Frequency, long long n, int d, vector<double> Prob) {
        double ChiSum = 0;
        for (int i = 0; i < d; ++i) {
            ChiSum += ((double)Frequency[i] - Prob[i] * n) * ((double)Frequency[i] - Prob[i] * n) / (Prob[i] * n);
        }
        return ChiSum;
    }

    int Comparator(int a, int b, int c, int d, int e) {
        vector <int> Counter(10, 0);
        vector <int> Values{ a,b,c,d,e };
        int Sum = 0;
        int i = 0;
        for (i = 0; i < 5; ++i) {
            if (Counter[Values[i]] == 0) {
                Counter[Values[i]]++;
                Sum++;
            }
        }
        return (Sum - 1);
    }
};

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

int main() {
    setlocale(LC_ALL, "Russian");
    NumberGenerator Generator = NumberGenerator(4294967296, 2147001325, 715136305, 0);
    Criteria Checker;
    double V, ChiQuant, alpha;
    cout << "1 -- проверка по критерию частот" << endl;
    cout << "2 -- проверка по критерию серий" << endl;
    cout << "3 -- проверка по критерию интервалов" << endl;
    cout << "4 -- проверка по покер-критерию" << endl;
    cout << "остальное -- выход" << endl;
    cout << "Введите номер команды: ";
    int CommandNumber;
    cin >> CommandNumber;
    switch (CommandNumber) {
        case 1:
            V = Checker.FrequencyTest(Generator, 100, 2147483648);
            cout << "Статистика V равна: " << V << endl;
            cout << "Введите альфу, для которой хотите проверить гипотезу: ";
            cin >> alpha;
            ChiQuant = GoldSteinChiSquare(alpha, 100-1);
            cout << "Хи-квадрат для альфа равное: " << ChiQuant << endl;
            if (V < ChiQuant) {
                cout << "Гипотеза принимается" << endl;
            }
            else {
                cout << "Гипотеза отвергается" << endl;
            }
            break;
        case 2:
            V = Checker.SerialTest(Generator, 100, 2147483648/2);
            cout << "Статистика V равна: " << V << endl;
            cout << "Введите альфу, для которой хотите проверить гипотезу: ";
            cin >> alpha;
            ChiQuant = GoldSteinChiSquare(alpha, (100*100)-1);
            cout << "Хи-квадрат для альфа равное: " << ChiQuant << endl;
            if (V < ChiQuant) {
                cout << "Гипотеза принимается" << endl;
            }
            else {
                cout << "Гипотеза отвергается" << endl;
            }
            break;
        case 3:
            V = Checker.CouponTest(Generator, 10, 2147483648);
            cout << "Статистика V равна: " << V << endl;cout << "Статистика V равна: " << V << endl;
            cout << "Введите альфу, для которой хотите проверить гипотезу: ";
            cin >> alpha;
            ChiQuant = GoldSteinChiSquare(alpha, 5);
            cout << "Хи-квадрат для альфа равное: " << ChiQuant << endl;
            if (V < ChiQuant) {
                cout << "Гипотеза принимается" << endl;
            }
            else {
                cout << "Гипотеза отвергается" << endl;
            }
            break;
        case 4:
            V = Checker.PokerTest(Generator, 10, 12333+2); // чтобы делилось на 5
            cout << "Статистика V равна: " << V << endl;
            cout << "Введите альфу, для которой хотите проверить гипотезу: ";
            cin >> alpha;
            ChiQuant = GoldSteinChiSquare(alpha, 4);
            cout << "Хи-квадрат для альфа равное: " << ChiQuant << endl;
            if (V < ChiQuant) {
                cout << "Гипотеза принимается" << endl;
            }
            else {
                cout << "Гипотеза отвергается" << endl;
            }
            break;
        default:
            return 0;
    }
}
