#include <iostream>
#include <cmath>
#include <fstream>

const double e = 2.7;

void iskhodsystem(double x, double v, double& dxdt, double& dvdt) {  //исходная системка
    dxdt = v; 
    dvdt = (-8 * x * pow(e, x * x)) / 9.1;
}

// Метод Эйлера
void eulerMethod(double x, double v, double h, double& xNew, double& vNew) {
    double dxdt, dvdt;
    iskhodsystem(x, v, dxdt, dvdt);  //производные
    xNew = x + h * dxdt; 
    vNew = v + h * dvdt;  
}


int main() {
    double tStart = 0;  //начальное t
    double tEnd = 15;   //конец t

    double x0 = 0;
    double v0 = 1;

    //Шаги интегрирования
    double hValues[] = { 1, 0.1, 0.01, 0.001 };
    const int nH = sizeof(hValues) / sizeof(hValues[0]);

    std::string fileNames[] = { "euler1.txt", "euler01.txt", "euler001.txt", "euler0001.txt" };

    for (int i = 0; i < nH; ++i) {
        double h = hValues[i];
        int nSteps = static_cast<int>((tEnd - tStart) / h);

        double* tValues = new double[nSteps + 1];
        double* xValues = new double[nSteps + 1];
        double* vValues = new double[nSteps + 1];

        // Начальные условия
        tValues[0] = tStart;
        xValues[0] = x0;
        vValues[0] = v0;

        for (int j = 0; j < nSteps; ++j) {
            double xNew, vNew;
            eulerMethod(xValues[j], vValues[j], h, xNew, vNew);
            tValues[j + 1] = tValues[j] + h;
            xValues[j + 1] = xNew;
            vValues[j + 1] = vNew;
        }

        std::ofstream outputFile(fileNames[i]);
        if (outputFile.is_open()) {
            for (int j = 0; j <= nSteps; ++j) {
                outputFile << tValues[j] << " " << xValues[j] << " " << vValues[j] << std::endl;
            }
            outputFile.close();
            std::cout << "Результаты для h = " << h << " записаны в файл " << fileNames[i] << std::endl;
        }
        else {
            std::cout << "Не удалось открыть файл " << fileNames[i] << " для записи" << std::endl;
        }

        delete[] tValues;
        delete[] xValues;
        delete[] vValues;
    }

    return 0;
}
