
#include <iostream>
#include <cmath>
#include <fstream>

//������� ���������������� ���������
void system(double state[2], double t, double dxdt[2]) {
    double x = state[0];
    double v = state[1];

    dxdt[0] = v; // dx/dt = v
    dxdt[1] = -8 / 9.1 * x * std::exp(x * x); // dv/dt = -8/9.1 * x * exp(x^2)
}

//������ �����-����� 4-�� �������
void rungeKutta4(double state[2], double t, double dt) {
    double k1[2], k2[2], k3[2], k4[2];

    // ���������� k1
    system(state, t, k1);

    // ���������� k2
    double tempState[2] = { state[0] + 0.5 * dt * k1[0], state[1] + 0.5 * dt * k1[1] };
    system(tempState, t + 0.5 * dt, k2);

    // ���������� k3
    tempState[0] = state[0] + 0.5 * dt * k2[0];
    tempState[1] = state[1] + 0.5 * dt * k2[1];
    system(tempState, t + 0.5 * dt, k3);

    // ���������� k4
    tempState[0] = state[0] + dt * k3[0];
    tempState[1] = state[1] + dt * k3[1];
    system(tempState, t + dt, k4);

    // ���������� ���������
    state[0] += (dt / 6.0) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
    state[1] += (dt / 6.0) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
}

int main() {

    setlocale(LC_ALL, "RUS");
    double initialState[2] = { 0.0, 1.0 }; // ��������� ��������� x=0, v=1
    double endTime = 15.0; // �������� ������ �������

    double dtValues[] = { 1.0, 0.1, 0.01, 0.001 };
    int numSteps = sizeof(dtValues) / sizeof(dtValues[0]);

    for (int i = 0; i < numSteps; ++i) {
        double t = 0.0; // ��������� ������ �������
        double dt = dtValues[i]; // ������� ��� �� �������

        // �������� ����� �����
        std::string filename = "RK";
        if (dt == 1.0) filename += "1";
        else if (dt == 0.1) filename += "01";
        else if (dt == 0.01) filename += "001";
        else if (dt == 0.001) filename += "0001";
        filename += ".txt";

        // �������� ����� ��� ������
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "������ ��� �������� ����� " << filename << std::endl;
            return 1;
        }

        // ����������� ���������� ��������� ��� ������� ����
        double currentState[2] = { initialState[0], initialState[1] };

        while (t < endTime) {
            file << t << " " << currentState[0] << " " << currentState[1] << std::endl;

            rungeKutta4(currentState, t, dt);
            t += dt;
        }

        file.close();
    }

    // �������� ������ SKIRK ��� ������� ����
    for (int i = 0; i < numSteps; ++i) {
        double dt = dtValues[i]; // ������� ��� �� �������

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
            std::cerr << "������ ��� �������� ������ �� ������." << std::endl;
            return 1;
        }

        // �������� ����� ��� ������
        std::ofstream fileSKIRK(filenameSKIRK);
        if (!fileSKIRK.is_open()) {
            std::cerr << "������ ��� �������� ����� ��� ������." << std::endl;
            return 1;
        }

        double t, x, v; // ��� ����� RK
        double skipT, skipX, skipV; // ��� ����� SKIPYDATA

        while (fileRK >> t >> x >> v && fileSKIP >> skipT >> skipX >> skipV) {
            // ���������� ������� �� ������ ��� ������� � �������� ��������
            double diffX = std::abs(x - skipX);
            double diffV = std::abs(v - skipV);

            // ������ ����������� � ����
            fileSKIRK << t << " " << diffX << " " << diffV << std::endl;
        }

        fileRK.close();
        fileSKIP.close();
        fileSKIRK.close();
    }

    return 0;
}


