#include "Main_Header.h"

// Векторная функция F
vector<double> F(const vector<double>& x)
{
    double f1 = tan(x[0] * x[1] + 0.2) - pow(x[0], 2);
    double f2 = 0.5 * pow(x[0], 2) + 2 * pow(x[1], 2) - 1;
    return { f1, f2 };
}

// Аналитическая матрица Якоби
vector<vector<double>> J_analytical(const vector<double>& x)
{
    double x1 = x[0], x2 = x[1];
    double sec_sq = 1.0 / pow(cos(x1 * x2 + 0.2), 2);
    return {
        {x2 * sec_sq - 2 * x1, x1 * sec_sq},
        {x1, 4 * x2}
    };
}

// Численная матрица Якоби методом конечных разностей
vector<vector<double>> J_numerical(const vector<double>& x, double M)
{
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

// Функция, реализующая метод Ньютона.
// Если параметр M < 0, используется аналитическое вычисление Якобиана, иначе — численное с указанным M.
void NewtonMethod(const vector<double>& x0, double eps1, double eps2, int max_iter, double M) {
    vector<double> x = x0;
    cout << "| Итерация | delta1 (||F(x)||) | delta2 (норма Δx) |" << endl;
    cout << "|----------|-------------------|-------------------|" << endl;

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

        cout << "| " << setw(8) << k + 1 << " | " << scientific << setprecision(11) << delta1 << " | " << setprecision(11) << delta2 << " |" << endl;

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