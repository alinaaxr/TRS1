#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

//Система
std::pair<double, double> system(double x, double v) {
    double dxdt = v;
    double dvdt = -8 / 9.1 * x * exp(x * x);
    return std::make_pair(dxdt, dvdt);
}

// Явная двухшаговая схема Адамса
std::pair<double, double> adamsBashforth2(double x_prev, double v_prev, double x_curr, double v_curr, double h) {
    double dxdt_curr = v_curr;  //производные
    double dvdt_curr = -8 / 9.1 * x_curr * exp(x_curr * x_curr);

    double dxdt_prev = v_prev;  // производные в предыдущей точке
    double dvdt_prev = -8 / 9.1 * x_prev * exp(x_prev * x_prev);

    // Применяем схему Адамса-Бэшфорта второго порядка
    double x_next = x_curr + h * (3 / 2.0 * dxdt_curr - 1 / 2.0 * dxdt_prev);
    double v_next = v_curr + h * (3 / 2.0 * dvdt_curr - 1 / 2.0 * dvdt_prev);

    return std::make_pair(x_next, v_next);
}

// Метод Рунге-Кутта 4-го порядка
std::pair<double, double> rungeKutta4(double x, double v, double t, double h) {
    double k1_x = h * v;
    double k1_v = h * (-8 / 9.1 * x * exp(x * x));

    double k2_x = h * (v + 0.5 * k1_v);
    double k2_v = h * (-8 / 9.1 * (x + 0.5 * k1_x) * exp((x + 0.5 * k1_x) * (x + 0.5 * k1_x)));

    double k3_x = h * (v + 0.5 * k2_v);
    double k3_v = h * (-8 / 9.1 * (x + 0.5 * k2_x) * exp((x + 0.5 * k2_x) * (x + 0.5 * k2_x)));

    double k4_x = h * (v + k3_v);
    double k4_v = h * (-8 / 9.1 * (x + k3_x) * exp((x + k3_x) * (x + k3_x)));

    double x_new = x + (1 / 6.0) * (k1_x + 2 * k2_x + 2 * k3_x + k4_x);
    double v_new = v + (1 / 6.0) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v);

    return std::make_pair(x_new, v_new);
}

int main() {
    // Параметры симуляции
    double t_start = 0;
    double t_end = 15;

    // Различные шаги интегрирования
    double h_values[] = { 1.0, 0.1, 0.01, 0.001 };
    std::string file_names[] = { "adams1.txt", "adams01.txt", "adams001.txt", "adams0001.txt" };

    for (int j = 0; j < 4; ++j) {
        double h = h_values[j];
        int n = static_cast<int>((t_end - t_start) / h) + 1;

        // Инициализация массивов для хранения результатов
        std::vector<double> t_values(n);
        std::vector<double> x_values(n);
        std::vector<double> v_values(n);

        // Начальные условия
        x_values[0] = 0;
        v_values[0] = 1;

        // Заполнение массива t_values
        for (int i = 0; i < n; ++i) {
            t_values[i] = t_start + i * h;
        }

        // Вычисляем вторую точку методом Рунге-Кутта
        auto new_state = rungeKutta4(x_values[0], v_values[0], t_values[0], h);
        x_values[1] = new_state.first;
        v_values[1] = new_state.second;

        // Симуляция
        for (int i = 1; i < n - 1; ++i) {
            auto new_state = adamsBashforth2(x_values[i - 1], v_values[i - 1], x_values[i], v_values[i], h);
            x_values[i + 1] = new_state.first;
            v_values[i + 1] = new_state.second;
        }

        // Обработка последней точки методом Рунге-Кутта
        if (n > 1) {
            auto new_state = rungeKutta4(x_values[n - 2], v_values[n - 2], t_values[n - 2], h);
            x_values[n - 1] = new_state.first;
            v_values[n - 1] = new_state.second;
        }

        // Визуализация результатов (записываем в файл для последующего построения графика)
        std::ofstream file(file_names[j]);
        if (file.is_open()) {
            for (int i = 0; i < n; ++i) {
                file << t_values[i] << " " << x_values[i] << " " << v_values[i] << "\n";
            }
            file.close();
        }
        else {
            std::cout << "Unable to open file";
        }
    }

    // Объявляем переменные и потоки один раз
    std::string adamsFile = "adams0001.txt";
    std::string rungeFile = "runge0001.txt";
    std::string eulerFile = "euler0001.txt";
    std::string outputRungeFile = "aDRU0001.txt";
    std::string outputFile = "AdEu0001.txt";

    std::vector<double> adamsT, adamsX, adamsV;
    std::vector<double> rungeX, rungeV;
    std::vector<double> eulerX, eulerV;

    // Чтение данных из файла Adams0001.txt
    std::ifstream adamsStream(adamsFile);
    if (adamsStream.is_open()) {
        double t, x, v;
        while (adamsStream >> t >> x >> v) {
            adamsT.push_back(t);
            adamsX.push_back(x);
            adamsV.push_back(v);
        }
        adamsStream.close();
    }
    else {
        std::cout << "Unable to open Adams file";
        return 1;
    }

    // Чтение данных из файла Runge0001.txt
    std::ifstream rungeStream(rungeFile);
    if (rungeStream.is_open()) {
        double t, x, v;
        while (rungeStream >> t >> x >> v) {
            rungeX.push_back(x);
            rungeV.push_back(v);
        }
        rungeStream.close();
    }
    else {
        std::cout << "Unable to open Runge file";
        return 1;
    }

    // Чтение данных из файла Euler0001.txt
    std::ifstream eulerStream(eulerFile);
    if (eulerStream.is_open()) {
        double t, x, v;
        while (eulerStream >> t >> x >> v) {
            eulerX.push_back(x);
            eulerV.push_back(v);
        }
        eulerStream.close();
    }
    else {
        std::cout << "Unable to open Euler file";
        return 1;
    }

    // Проверка размеров векторов
    if (adamsX.size() != rungeX.size() || adamsV.size() != rungeV.size() ||
        adamsX.size() != eulerX.size() || adamsV.size() != eulerV.size()) {
        std::cout << "Sizes of files do not match";
        return 1;
    }

    // Запись результатов в новый файл для Adams и Runge
    std::ofstream outputRunge(outputRungeFile);
    if (outputRunge.is_open()) {
        for (size_t i = 0; i < adamsT.size(); ++i) {
            double diffX = std::abs(adamsX[i] - rungeX[i]);
            double diffV = std::abs(adamsV[i] - rungeV[i]);
            outputRunge << adamsT[i] << " " << diffX << " " << diffV << "\n";
        }
        outputRunge.close();
    }
    else {
        std::cout << "Unable to open output file";
        return 1;
    }

    // Запись результатов в новый файл для Adams и Euler
    std::ofstream output(outputFile);
    if (output.is_open()) {
        for (size_t i = 0; i < adamsT.size(); ++i) {
            double diffX = std::abs(adamsX[i] - eulerX[i]);
            double diffV = std::abs(adamsV[i] - eulerV[i]);
            output << adamsT[i] << " " << diffX << " " << diffV << "\n";
        }
        output.close();
    }
    else {
        std::cout << "Unable to open output file";
        return 1;
    }

    return 0;
}
