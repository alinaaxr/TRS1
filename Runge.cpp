//#include <iostream>
//#include <cmath>
//#include <fstream>
//#include <vector>
//
//// Функция, представляющая систему дифференциальных уравнений
//void system(double x, double v, double& dxdt, double& dvdt) {
//    dxdt = v;
//    dvdt = -8 / 9.1 * x * std::exp(x * x);
//}
//
//// Метод Рунге-Кутта 4-го порядка
//void rungeKutta4(double x, double v, double h, double& xNew, double& vNew) {
//    double k1_x, k1_v, k2_x, k2_v, k3_x, k3_v, k4_x, k4_v;
//
//    // Вычисление k1
//    double dxdt, dvdt;
//    system(x, v, dxdt, dvdt);
//    k1_x = h * dxdt;
//    k1_v = h * dvdt;
//
//    // Вычисление k2
//    system(x + 0.5 * k1_x, v + 0.5 * k1_v, dxdt, dvdt);
//    k2_x = h * dxdt;
//    k2_v = h * dvdt;
//
//    // Вычисление k3
//    system(x + 0.5 * k2_x, v + 0.5 * k2_v, dxdt, dvdt);
//    k3_x = h * dxdt;
//    k3_v = h * dvdt;
//
//    // Вычисление k4
//    system(x + k3_x, v + k3_v, dxdt, dvdt);
//    k4_x = h * dxdt;
//    k4_v = h * dvdt;
//
//    // Вычисление новых значений x и v
//    xNew = x + (1.0 / 6.0) * (k1_x + 2 * k2_x + 2 * k3_x + k4_x);
//    vNew = v + (1.0 / 6.0) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v);
//}
//
//// Функция для чтения данных из файла и вычисления разницы
//void calculateDifference(const std::string& eulerFile, const std::string& rungeFile, const std::string& output) {
//    std::ifstream eulerInput(eulerFile);
//    std::ifstream rungeInput(rungeFile);
//    std::ofstream outputFile(output);
//
//    if (eulerInput.is_open() && rungeInput.is_open() && outputFile.is_open()) {
//        double eulerT, eulerX, eulerV;
//        double rungeT, rungeX, rungeV;
//
//        while (eulerInput >> eulerT >> eulerX >> eulerV && rungeInput >> rungeT >> rungeX >> rungeV) {
//            double differenceX = std::abs(eulerX - rungeX);
//            double differenceV = std::abs(eulerV - rungeV);
//            outputFile << eulerT << " " << differenceX << " " << differenceV << std::endl;
//        }
//
//        eulerInput.close();
//        rungeInput.close();
//        outputFile.close();
//
//        std::cout << "Разница записана в файл " << output << std::endl;
//    }
//    else {
//        std::cout << "Не удалось открыть один из файлов для чтения или записи." << std::endl;
//    }
//}
//
//
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
//    std::string fileNames[] = { "Runge1.txt", "Runge01.txt", "Runge001.txt", "Runge0001.txt" };
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
//            rungeKutta4(xValues[j], vValues[j], h, xNew, vNew);
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
//    // Вычисление разницы между файлами
//    // Вычисление разницы между файлами
//    calculateDifference("EULER1.txt", "Runge1.txt", "RAZN1.txt");
//    calculateDifference("EULER01.txt", "Runge01.txt", "RAZN01.txt");
//    calculateDifference("EULER001.txt", "Runge001.txt", "RAZN001.txt");
//    calculateDifference("EULER0001.txt", "Runge0001.txt", "RAZN0001.txt");
//
//
//    return 0;
//}
#include <iostream>
#include <cmath>
#include <fstream>

// Функция, описывающая систему дифференциальных уравнений
void system(double state[2], double t, double dxdt[2]) {
    double x = state[0];
    double v = state[1];

    dxdt[0] = v; // dx/dt = v
    dxdt[1] = -8 / 9.1 * x * std::exp(x * x); // dv/dt = -8/9.1 * x * exp(x^2)
}

// Функция метода Рунге-Кутты 4-го порядка
void rungeKutta4(double state[2], double t, double dt) {
    double k1[2], k2[2], k3[2], k4[2];

    // Вычисление k1
    system(state, t, k1);

    // Вычисление k2
    double tempState[2] = { state[0] + 0.5 * dt * k1[0], state[1] + 0.5 * dt * k1[1] };
    system(tempState, t + 0.5 * dt, k2);

    // Вычисление k3
    tempState[0] = state[0] + 0.5 * dt * k2[0];
    tempState[1] = state[1] + 0.5 * dt * k2[1];
    system(tempState, t + 0.5 * dt, k3);

    // Вычисление k4
    tempState[0] = state[0] + dt * k3[0];
    tempState[1] = state[1] + dt * k3[1];
    system(tempState, t + dt, k4);

    // Обновление состояния
    state[0] += (dt / 6.0) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
    state[1] += (dt / 6.0) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
}

int main() {

    setlocale(LC_ALL, "RUS");
    double initialState[2] = { 0.0, 1.0 }; // Начальное состояние x=0, v=1
    double endTime = 15.0; // Конечный момент времени

    // Шаги по времени
    double dtValues[] = { 1.0, 0.1, 0.01, 0.001 };
    int numSteps = sizeof(dtValues) / sizeof(dtValues[0]);

    for (int i = 0; i < numSteps; ++i) {
        double t = 0.0; // Начальный момент времени
        double dt = dtValues[i]; // Текущий шаг по времени

        // Создание имени файла
        std::string filename = "RK";
        if (dt == 1.0) filename += "1";
        else if (dt == 0.1) filename += "01";
        else if (dt == 0.01) filename += "001";
        else if (dt == 0.001) filename += "0001";
        filename += ".txt";

        // Открытие файла для записи
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Ошибка при открытии файла " << filename << std::endl;
            return 1;
        }

        // Копирование начального состояния для каждого шага
        double currentState[2] = { initialState[0], initialState[1] };

        while (t < endTime) {
            file << t << " " << currentState[0] << " " << currentState[1] << std::endl;

            rungeKutta4(currentState, t, dt);
            t += dt;
        }

        file.close();
    }

    // Создание файлов SKIRK для каждого шага
    for (int i = 0; i < numSteps; ++i) {
        double dt = dtValues[i]; // Текущий шаг по времени

        std::string filenameRK = "RK";
        if (dt == 1.0) filenameRK += "1";
        else if (dt == 0.1) filenameRK += "01";
        else if (dt == 0.01) filenameRK += "001";
        else if (dt == 0.001) filenameRK += "0001";
        filenameRK += ".txt";

        std::string filenameSKIRK = "SKIRK";
        if (dt == 1.0) filenameSKIRK += "1";
        else if (dt == 0.1) filenameSKIRK += "01";
        else if (dt == 0.01) filenameSKIRK += "001";
        else if (dt == 0.001) filenameSKIRK += "0001";
        filenameSKIRK += ".txt";

        std::ifstream fileRK(filenameRK);
        std::ifstream fileSKIP("SKIPYDATA001.txt");
        if (!fileRK.is_open() || !fileSKIP.is_open()) {
            std::cerr << "Ошибка при открытии одного из файлов." << std::endl;
            return 1;
        }

        // Открытие файла для записи
        std::ofstream fileSKIRK(filenameSKIRK);
        if (!fileSKIRK.is_open()) {
            std::cerr << "Ошибка при открытии файла для записи." << std::endl;
            return 1;
        }

        double t, x, v; // Для файла RK
        double skipT, skipX, skipV; // Для файла SKIPYDATA

        while (fileRK >> t >> x >> v && fileSKIP >> skipT >> skipX >> skipV) {
            // Вычисление разницы по модулю для второго и третьего столбцов
            double diffX = std::abs(x - skipX);
            double diffV = std::abs(v - skipV);

            // Запись результатов в файл
            fileSKIRK << t << " " << diffX << " " << diffV << std::endl;
        }

        fileRK.close();
        fileSKIP.close();
        fileSKIRK.close();
    }

    return 0;
}


