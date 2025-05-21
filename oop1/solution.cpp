#include "Solution.h"

Solution::Solution(double A, double B, double C) : a(A), b(B), c(C) {}

void Solution::print() {
    cout << "Коэффициенты квадратного уравнения: a = " << a << ", b = " << b << ", c = " << c << endl;
}

void Solution::solve() {
    double D = b * b - 4 * a * c;

    if (D > 0) {
        double x1 = (-b + sqrt(D)) / (2 * a);
        double x2 = (-b - sqrt(D)) / (2 * a);
        cout << "Корни уравнения: x1 = " << x1 << ", x2 = " << x2 << endl;
    } else if (D == 0) {
        double x = -b / (2 * a);
        cout << "Единственный корень: x = " << x << endl;
    } else {
        cout << "Корней нет (D < 0, комплексные числа)." << endl;
    }
}

string Solution::rootType() {
    double D = b * b - 4 * a * c;
    if (D > 0) return "Два вещественных корня";
    else if (D == 0) return "Один вещественный корень";
    else return "Комплексные корни";
}

void Solution::readFromFile(const string& filename) {
    ifstream file(filename);
    if (file.is_open()) {
        file >> a >> b >> c;
        file.close();
    } else {
        cout << "Ошибка: не удалось открыть файл!" << endl;
    }
}
