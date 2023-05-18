#include <cfloat>
#include <fstream>
#include <iostream>
#include <math.h>
#include <time.h>
#include <vector>
#include <sstream>
#include <chrono>


struct S_Monkey {
    std::vector<double> c;  // coordinates
    std::vector<double> cB; // best coordinates
    double h;               // height of the mountain
    double hB;              // best height of the mountain
    int lCNT;               // local search counter
};

class C_AO_MA {
    //----------------------------------------------------------------------------
public:
    std::vector<S_Monkey> m;       // monkeys
    std::vector<double> rangeMax;  // maximum search range
    std::vector<double> rangeMin;  // manimum search range
    std::vector<double> rangeStep; // step search
    std::vector<double> cB;        // best coordinates
    double hB;                     // best height of the mountain

    void Init(const int coordNumberP, // coordinates number
        const int monkeysNumberP,     // monkeys number
        const double bCoefficientP,   // local search coefficient
        const double vCoefficientP,   // jump coefficient
        const int jumpsNumberP);      // jumps number

    void Moving();
    void Revision();

    //----------------------------------------------------------------------------
private:
    int coordNumber;   // coordinates number
    int monkeysNumber; // monkeys number

    std::vector<double> b;  // local search coefficient
    std::vector<double> v;  // jump coefficient
    double bCoefficient;    // local search coefficient
    double vCoefficient;    // jump coefficient
    double jumpsNumber;     // jumps number
    std::vector<double> cc; // coordinate center

    bool revision;

    double SeInDiSp(double In, double InMin, double InMax, double Step);
    double RNDfromCI(double min, double max);
};

void C_AO_MA::Init(const int coordNumberP, // coordinates number
    const int monkeysNumberP,              // monkeys number
    const double bCoefficientP,            // local search coefficient
    const double vCoefficientP,            // jump coefficient
    const int jumpsNumberP)                // jumps number
{
    // MathSrand((int)GetMicrosecondCount()); // reset of the generator
    hB = -DBL_MAX;
    revision = false;

    coordNumber = coordNumberP;
    monkeysNumber = monkeysNumberP;
    bCoefficient = bCoefficientP;
    vCoefficient = vCoefficientP;
    jumpsNumber = jumpsNumberP;

    rangeMax.resize(coordNumber);  //
    rangeMin.resize(coordNumber);  //
    rangeStep.resize(coordNumber); //
    b.resize(coordNumber);         //
    v.resize(coordNumber);         //
    cc.resize(coordNumber);        //

    m.resize(monkeysNumber);

    for (int i = 0; i < monkeysNumber; i++) {
        m[i].c.resize(coordNumber);  //
        m[i].cB.resize(coordNumber); //
        m[i].h = -DBL_MAX;
        m[i].hB = -DBL_MAX;
        m[i].lCNT = 0;
    }

    cB.resize(coordNumber); //
}

void C_AO_MA::Moving()
{
    //----------------------------------------------------------------------------
    if (!revision) {
        hB = -DBL_MAX;

        for (int monk = 0; monk < monkeysNumber; monk++) {
            for (int c = 0; c < coordNumber; c++) {
                m[monk].c[c] = RNDfromCI(rangeMin[c], rangeMax[c]);
                m[monk].c[c] = SeInDiSp(m[monk].c[c], rangeMin[c], rangeMax[c], rangeStep[c]);
                m[monk].h = -DBL_MAX;
                m[monk].hB = -DBL_MAX;
                m[monk].lCNT = 0;
            }
        }

        for (int c = 0; c < coordNumber; c++) {
            v[c] = (rangeMax[c] - rangeMin[c]) * vCoefficient;
            b[c] = (rangeMax[c] - rangeMin[c]) * bCoefficient;
        }

        revision = true;
    }
    //----------------------------------------------------------------------------
    else {
        double r1 = 0.0;
        double r2 = 0.0;

        // calculate the coordinate center of the monkeys----------------------------
        for (int c = 0; c < coordNumber; c++) {
            cc[c] = 0.0;

            for (int monk = 0; monk < monkeysNumber; monk++) {
                cc[c] += m[monk].cB[c];
            }

            cc[c] /= monkeysNumber;
        }

        //==========================================================================
        // monkeys jumping ----------------------------------------------------------
        for (int monk = 0; monk < monkeysNumber; monk++) {
            if (m[monk].lCNT < jumpsNumber)
                // local jump--------------------------------------------------------------
            {
                for (int c = 0; c < coordNumber; c++) {
                    m[monk].c[c] = RNDfromCI(m[monk].cB[c] - b[c], m[monk].cB[c] + b[c]);
                    m[monk].c[c] = SeInDiSp(m[monk].c[c], rangeMin[c], rangeMax[c], rangeStep[c]);
                }
            }
            else {
                // global jump-------------------------------------------------------------
                for (int c = 0; c < coordNumber; c++) {
                    r1 = RNDfromCI(0.0, 1.0);
                    r1 = r1 > 0.5 ? 1.0 : -1.0;
                    r2 = RNDfromCI(1.0, 20.0);

                    // m [monk].c [c] = cc [c] + r1 * v [c] * pow (r2, -2.0);
                    m[monk].c[c] = cc[c] + v[c] * pow(r2, -2.0);

                    if (m[monk].c[c] < rangeMin[c])
                        m[monk].c[c] = rangeMax[c] - (rangeMin[c] - m[monk].c[c]);
                    if (m[monk].c[c] > rangeMax[c])
                        m[monk].c[c] = rangeMin[c] + (m[monk].c[c] - rangeMax[c]);

                    m[monk].c[c] = SeInDiSp(m[monk].c[c], rangeMin[c], rangeMax[c], rangeStep[c]);
                }
            }

            m[monk].lCNT++;
        }
    }
}

void C_AO_MA::Revision()
{
    for (int monk = 0; monk < monkeysNumber; monk++) {
        if (m[monk].h > hB) {
            hB = m[monk].h;
            // ArrayCopy(cB, m[monk].c, 0, 0, WHOLE_ARRAY);
            std::copy(m[monk].c.begin(), m[monk].c.end(), std::back_inserter(cB));
        }

        if (m[monk].lCNT <= jumpsNumber) // local jump
        {
            if (m[monk].h > m[monk].hB) {
                m[monk].hB = m[monk].h;
                // ArrayCopy(m[monk].cB, m[monk].c, 0, 0, WHOLE_ARRAY);
                std::copy(m[monk].c.begin(), m[monk].c.end(), std::back_inserter(m[monk].cB));
                m[monk].lCNT = 0;
            }
        }
        else {
            m[monk].hB = m[monk].h;
            // ArrayCopy(m[monk].cB, m[monk].c, 0, 0, WHOLE_ARRAY);
            std::copy(m[monk].c.begin(), m[monk].c.end(), std::back_inserter(m[monk].cB));
            m[monk].lCNT = 0;
        }
    }
}

//  Choice in discrete space
double C_AO_MA::SeInDiSp(double In, double InMin, double InMax, double Step)
{
    if (In <= InMin)
        return (InMin);
    if (In >= InMax)
        return (InMax);
    if (Step == 0.0)
        return (In);
    else
        return (InMin + Step * (double)round((In - InMin) / Step));
}

//  Random number generator in the custom interval.
double C_AO_MA::RNDfromCI(double min, double max)
{
    if (min == max)
        return (min);
    double Min, Max;
    if (min > max) {
        Min = max;
        Max = min;
    }
    else {
        Min = min;
        Max = max;
    }
    return (double(Min + ((Max - Min) * (double)rand() / 32767.0)));
}

// Scaling a number from a range to a specified range
double Scale(double In, double InMIN, double InMAX, double OutMIN, double OutMAX, bool Revers = false)
{
    if (OutMIN == OutMAX)
        return (OutMIN);
    if (InMIN == InMAX)
        return ((OutMIN + OutMAX) / 2.0);
    else {
        if (Revers) {
            if (In < InMIN)
                return (OutMAX);
            if (In > InMAX)
                return (OutMIN);
            return (((InMAX - In) * (OutMAX - OutMIN) / (InMAX - InMIN)) + OutMIN);
        }
        else {
            if (In < InMIN)
                return (OutMIN);
            if (In > InMAX)
                return (OutMAX);
            return (((In - InMIN) * (OutMAX - OutMIN) / (InMAX - InMIN)) + OutMIN);
        }
    }
}
std::string doubleToString(double number, int precision) {
    std::ostringstream stream;
    stream.precision(precision);
    stream << std::fixed << number;
    return stream.str();
}

std::vector<int> Population_P = {20, 50, 70, 100, 150};        // Population size
std::vector<double> Bcoefficient_P = {0.01, 0.03, 0.05, 0.07, 0.09}; // Local search coefficient
std::vector<double> Vcoefficient_P = {0.5, 0.6, 0.7, 0.8, 0.9};  // Jump coefficient
std::vector<int> JumpsNumber_P = {50, 60, 70, 80, 90};       // Jumps number

double ArgumentStep_P = 0.0;    // Argument Step
std::vector<int> NumbTestFuncRuns_P = { 1000, 10000, 100000}; // Number of test function runs
int NumberRepetTest_P = 5;      // Test repets number

C_AO_MA AO; // AO object

class function_proto {
public:
    virtual double func(std::vector<double>& v) = 0;
    virtual double GetMaxX(int i) = 0;
    virtual double GetMinX(int i) = 0;
    virtual double GetMaxFun() = 0;
    virtual double GetMinFun() = 0;
    virtual int count_x() = 0;
    virtual std::string get_name() = 0;
};

class eggholder : public function_proto {
public:
    eggholder(int sz, std::string nm) : name(nm)
    {
        Ccount_x = sz;
        MaxX.resize(sz);
        MinX.resize(sz);
        for (int i = 0; i < sz; i++) {
            MaxX[i] = 512;
            MinX[i] = -512;
        }
        MaxFun = 862.1047;
        MinFun = Mins[sz-2];
    }
    double func(std::vector<double>& v)
    {
        double res = 0;
        for (int i = 0; i < Ccount_x - 1; i++) {
           res -= (v[i+1] + 47) * std::sin(std::sqrt(std::abs(v[i+1] + 47 + (v[i] / 2.0)))) +
                v[i]*std::sin(std::sqrt(std::abs(v[i]-(v[i+1] + 47))));
        }
        return res;
    }
    double Mins[9]{ -959.6406, -1888.3213909, -2808.1847922, -3719.7248363, -4625.1447737,
     -5548.9775483, -6467.0193267, -7376.2797668, -8291.2400675};
    double GetMaxX(int i) { return MaxX[i]; }
    double GetMinX(int i) { return MinX[i]; }
    double GetMaxFun() { return MaxFun; }
    double GetMinFun() { return MinFun; }
    int count_x() { return Ccount_x; }
    std::string get_name() { return name; }

private:
    std::vector<double> MaxX;
    std::vector<double> MinX;
    double MaxFun;
    double MinFun;
    int Ccount_x;
    std::string name;
};

void FuncTests(function_proto& f, std::string fileName)
{
    std::ofstream out;
    out.open(fileName + ".txt");
    int xConv = 0.0;
    int yConv = 0.0;
    double aveResult = 0.0;
    for (int count_r = 0; count_r < NumbTestFuncRuns_P.size(); count_r++) {
        for (int PP = 0; PP < Population_P.size(); PP++) {
            for (int BP = 0; BP < Bcoefficient_P.size(); BP++) {
                for (int VP = 0; VP < Vcoefficient_P.size(); VP++) {
                    for (int JP = 0; JP < JumpsNumber_P.size(); JP++) {
                       //----------------------------------------------------------------------------
                        for (int test = 0; test < NumberRepetTest_P; test++) {
                            //--------------------------------------------------------------------------
                            int epochCount = NumbTestFuncRuns_P[count_r] / Population_P[PP];
                            AO.Init(f.count_x(), Population_P[PP], Bcoefficient_P[BP], Vcoefficient_P[VP], JumpsNumber_P[JP]);

                            for (int i = 0; i < f.count_x(); i++) {
                                AO.rangeMax[i] = f.GetMaxX(i);
                                AO.rangeMin[i] = f.GetMinX(i);
                                AO.rangeStep[i] = ArgumentStep_P;
                            }

                            // Optimization-------------------------------------------------------------
                            for (int epochCNT = 1; epochCNT <= epochCount; epochCNT++) {
                                AO.Moving();

                                for (int set = 0; set < AO.m.size(); set++) {
                                    AO.m[set].h = f.func(AO.m[set].c);
                                }

                                AO.Revision();
                            }

                            aveResult += AO.hB;
                         }

                        aveResult /= (double)NumberRepetTest_P;

                        double score = Scale(aveResult, f.GetMinFun(), f.GetMaxFun(), 0.0, 1.0, false);

                        if (out.is_open())
                        {
                            out << f.get_name() << "'s; Func runs " << NumbTestFuncRuns_P[count_r] 
                                << " Population: " << Population_P[PP] 
                                << " Local search coefficient: " << Bcoefficient_P[BP] 
                                << " Jump coefficient: " << Vcoefficient_P[VP] 
                                << " Jumps number: " << JumpsNumber_P[JP] 
                                << " result: " << -aveResult << std::endl;
                            out << "Score: " << doubleToString(score, 5) << std::endl;
                        }
                    }
                }
            }
        }
    }
    out.close();
}

int main(int argc, char* argv[])
{
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<double> a;
    a.resize(4);
    srand(time(NULL));
    for (int i = 2; i <= 10; i++) {
        eggholder f(i, "eggholder" + std::to_string(i));
        FuncTests(f, "eggholder" + std::to_string(i));
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Program execution time: " << duration << "seconds" << std::endl;
    return 0;
}
