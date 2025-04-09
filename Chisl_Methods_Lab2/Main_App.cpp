#include "Main_Header.h"

int main() {
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);

    vector<double> x0 = { 1.0, 1.0 };
    double eps1 = 1e-9, eps2 = 1e-9;
    int max_iter = 100;

    // Решение методом Ньютона с аналитическим Якобианом
    cout << "Аналитический метод:" << endl;
    NewtonMethod(x0, eps1, eps2, max_iter);

    // Решение методом Ньютона с численным вычислением Якобиана для разных значений M
    vector<double> M_values = { 0.01, 0.05, 0.1 };
    for (double M : M_values) {
        cout << "\nЧисленный метод (M = "<< setprecision(1) << M << "):" << endl;
        NewtonMethod(x0, eps1, eps2, max_iter, M);
    }

    return 0;
}
