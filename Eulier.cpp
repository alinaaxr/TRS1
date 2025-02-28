//#include <iostream>
//#include <cmath>
//#include <fstream>
//
//// Функция, представляющая систему дифференциальных уравнений
//void system(double x, double v, double& dxdt, double& dvdt) {
//    dxdt = v;
//    dvdt = -1 / 1.67 * (2 * x + 2 * x * std::pow(std::tan(x * x), 2));
//}
//
//// Метод Эйлера
//void eulerMethod(double x, double v, double h, double& xNew, double& vNew) {
//    double dxdt, dvdt;
//    system(x, v, dxdt, dvdt);
//    xNew = x + h * dxdt;
//    vNew = v + h * dvdt;
//}
//
//int main() {
//    // Параметры симуляции
//    double tStart = 0;
//    double tEnd = 15;
//
//    // Начальные условия
//    double x0 = 0;
//    double v0 = 1;
//
//    // Массив значений шага интегрирования
//    double hValues[] = { 1, 0.1, 0.01, 0.001 };
//    const int nH = sizeof(hValues) / sizeof(hValues[0]);
//
//    // Массив имён файлов для записи результатов
//    std::string fileNames[] = { "euler1.txt", "euler01.txt", "euler001.txt", "euler0001.txt" };
//
//    for (int i = 0; i < nH; ++i) {
//        double h = hValues[i];
//        int nSteps = static_cast<int>((tEnd - tStart) / h);
//
//        // Инициализация массивов для хранения результатов
//        double* tValues = new double[nSteps + 1];
//        double* xValues = new double[nSteps + 1];
//        double* vValues = new double[nSteps + 1];
//
//        // Начальные условия
//        tValues[0] = tStart;
//        xValues[0] = x0;
//        vValues[0] = v0;
//
//        // Симуляция
//        for (int j = 0; j < nSteps; ++j) {
//            double xNew, vNew;
//            eulerMethod(xValues[j], vValues[j], h, xNew, vNew);
//            tValues[j + 1] = tValues[j] + h;
//            xValues[j + 1] = xNew;
//            vValues[j + 1] = vNew;
//        }
//
//        // Визуализация результатов (вывод в файл для последующего построения графика)
//        std::ofstream outputFile(fileNames[i]);
//        if (outputFile.is_open()) {
//            for (int j = 0; j <= nSteps; ++j) {
//                outputFile << tValues[j] << " " << xValues[j] << " " << vValues[j] << std::endl;
//            }
//            outputFile.close();
//            std::cout << "Результаты для h = " << h << " записаны в файл " << fileNames[i] << std::endl;
//        }
//        else {
//            std::cout << "Не удалось открыть файл " << fileNames[i] << " для записи" << std::endl;
//        }
//
//        // Освобождение памяти
//        delete[] tValues;
//        delete[] xValues;
//        delete[] vValues;
//    }
//
//    return 0;
//}
