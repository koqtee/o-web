#include "Vector2.h"

int main() {
    // Создание объектов
    Vector2 v1(1, 2, 3, 4);
    Vector2 v2(4, 3, 2, 1);

    cout << "Вектор 1: " << v1 << endl;
    cout << "Вектор 2: " << v2 << endl;

    // Сложение векторов
    Vector2 vSum = v1 + v2;
    cout << "Сумма векторов: " << vSum << endl;

    // Умножение на число
    Vector2 vMult = v1 * 3;
    cout << "Вектор 1 * 3: " << vMult << endl;

    // Проверка оператора сравнения
    cout << "Векторы равны? " << (v1 == v2 ? "Да" : "Нет") << endl;

    // Унарный оператор минус
    Vector2 vNeg = -v1;
    cout << "Обратный знак Вектора 1: " << vNeg << endl;

    return 0;
}
