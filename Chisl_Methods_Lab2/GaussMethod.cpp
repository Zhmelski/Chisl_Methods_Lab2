#include "Main_Header.h"

// Решение СЛАУ методом Гаусса для системы 2×2
vector<double> GaussMethod(vector<vector<double>> A, vector<double> b) {
    int n = A.size();

    // Прямой ход (прямой эллиминация)
    for (int i = 0; i < n; ++i) {
        // Поиск максимального по модулю элемента в столбце i для устойчивости
        int p = i;
        for (int k = i + 1; k < n; ++k)
            if (fabs(A[k][i]) > fabs(A[p][i]))
                p = k;

        if (p != i)
        {
            swap(A[i], A[p]);
            swap(b[i], b[p]);
        }

        if (fabs(A[i][i]) < EPS)
            throw runtime_error("Нулевой элемент на главной диагонали");

        // Нормализация строки и исключение переменных из последующих строк
        for (int k = i + 1; k < n; ++k) {
            double factor = A[k][i] / A[i][i];
            b[k] -= factor * b[i];
            for (int j = i; j < n; ++j)
                A[k][j] -= factor * A[i][j];
        }
    }

    // Обратный ход (подстановка)
    vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[i];
        for (int j = i + 1; j < n; ++j)
            x[i] -= A[i][j] * x[j];
        x[i] /= A[i][i];
    }
    return x;
}