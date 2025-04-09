#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include <Windows.h>

using namespace std;

const double EPS = 1e-12;

// Векторная функция F
vector<double> F(const vector<double>& x) {
    double f1 = tan(x[0] * x[1] + 0.2) - pow(x[0], 2);
    double f2 = 0.5 * pow(x[0], 2) + 2 * pow(x[1], 2) - 1;
    return { f1, f2 };
}

// Аналитическая матрица Якоби
vector<vector<double>> J_analytical(const vector<double>& x) {
    double x1 = x[0], x2 = x[1];
    double sec_sq = 1.0 / pow(cos(x1 * x2 + 0.2), 2);
    return {
        {x2 * sec_sq - 2 * x1, x1 * sec_sq},
        {x1, 4 * x2}
    };
}

// Численная матрица Якоби методом конечных разностей
vector<vector<double>> J_numerical(const vector<double>& x, double M) {
    int n = x.size();
    vector<vector<double>> J(n, vector<double>(n, 0.0));
    vector<double> f = F(x);

    for (int j = 0; j < n; ++j) {
        vector<double> x_perturbed = x;
        double h = M * fabs(x[j]);
        h = (h < EPS) ? EPS : h; // Защита от деления на ноль
        x_perturbed[j] += h;

        vector<double> f_perturbed = F(x_perturbed);
        for (int i = 0; i < n; ++i) {
            J[i][j] = (f_perturbed[i] - f[i]) / h;
        }
    }
    return J;
}

// Решение СЛАУ методом Гаусса для системы 2×2
vector<double> GaussMethod(vector<vector<double>> A, vector<double> b) {
    int n = A.size();

    // Прямой ход
    for (int i = 0; i < n; ++i) {
        // Поиск максимального по модулю элемента в столбце i
        int p = i;
        for (int k = i + 1; k < n; ++k)
            if (fabs(A[k][i]) > fabs(A[p][i]))
                p = k;

        swap(A[i], A[p]);
        swap(b[i], b[p]);

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
        
    // Обратный ход
    vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[i];
        for (int j = i + 1; j < n; ++j)
            x[i] -= A[i][j] * x[j];
        x[i] /= A[i][i];
    }
    return x;
}

// Функция, реализующая метод Ньютона.
// Если параметр M < 0, используется аналитическое вычисление Якобиана, иначе — численное с указанным M.
void NewtonMethod(const vector<double>& x0, double eps1, double eps2, int max_iter, double M = -1.0) {
    vector<double> x = x0;
    cout << "| Итерация | beta1 (||F(x)||) | beta2 (норма Δx) |" << endl;
    cout << "|----------|------------------|------------------|" << endl;

    for (int k = 0; k < max_iter; ++k) {
        vector<double> f = F(x);
        double delta1 = max(fabs(f[0]), fabs(f[1]));

        vector<vector<double>> J;
        if (M < 0)
            J = J_analytical(x);
        else
            J = J_numerical(x, M);

        vector<double> delta_x;
        try {
            // Решаем систему J * delta_x = -F(x)
            delta_x = GaussMethod(J, vector<double>{ -f[0], -f[1] });
        }
        catch (const exception& e) {
            cerr << "Ошибка решения СЛАУ на итерации " << k + 1 << ": " << e.what() << endl;
            return;
        }

        // Новое приближение
        vector<double> x_new = { x[0] + delta_x[0], x[1] + delta_x[1] };

        // Вычисление δ₂ по относительной разности
        double delta2 = 0.0;
        for (int i = 0; i < 2; ++i) {
            double denom = (fabs(x_new[i]) >= 1.0) ? fabs(x_new[i]) : 1.0;
            delta2 = max(delta2, fabs((x_new[i] - x[i]) / denom));
        }

        cout << "| " << setw(8) << k + 1 << " | " << scientific << setprecision(10) << delta1 << " | " << setprecision(10) << delta2 << " |" << endl;

        // Проверка условий останова
        if (delta1 <= eps1 && delta2 <= eps2) {
            cout << "\nУспешная сходимость за " << k + 1 << " итераций." << endl;
            cout << "Найденное решение: x = {" << x_new[0] << ", " << x_new[1] << "}" << endl;
            return;
        }
        x = x_new;
    }
    cout << "\nДостигнуто максимальное число итераций." << endl;
    cout << "Последнее приближение: x = {" << x[0] << ", " << x[1] << "}" << endl;
}

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
        cout << "\nЧисленный метод (M = " << M << "):" << endl;
        NewtonMethod(x0, eps1, eps2, max_iter, M);
    }

    return 0;
}
