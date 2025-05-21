#include "Vector2.h"

// Конструктор по умолчанию
Vector2::Vector2() {
    for (int i = 0; i < 4; i++) {
        elements[i] = 0;
    }
}

// Конструктор с параметрами
Vector2::Vector2(int e1, int e2, int e3, int e4) {
    elements[0] = e1;
    elements[1] = e2;
    elements[2] = e3;
    elements[3] = e4;
}

// Метод печати
void Vector2::print() const {
    cout << "[" << elements[0] << ", " << elements[1] << ", " << elements[2] << ", " << elements[3] << "]" << endl;
}

// Перегрузка оператора сложения (+)
Vector2 Vector2::operator+(const Vector2& other) const {
    return Vector2(
        elements[0] + other.elements[0],
        elements[1] + other.elements[1],
        elements[2] + other.elements[2],
        elements[3] + other.elements[3]
    );
}

// Перегрузка оператора умножения на число (*)
Vector2 Vector2::operator*(int scalar) const {
    return Vector2(
        elements[0] * scalar,
        elements[1] * scalar,
        elements[2] * scalar,
        elements[3] * scalar
    );
}

// Перегрузка оператора сравнения (==)
bool Vector2::operator==(const Vector2& other) const {
    for (int i = 0; i < 4; i++) {
        if (elements[i] != other.elements[i]) {
            return false;
        }
    }
    return true;
}

// Перегрузка унарного оператора минус (-)
Vector2 Vector2::operator-() const {
    return Vector2(-elements[0], -elements[1], -elements[2], -elements[3]);
}

// Перегрузка оператора вывода (<<)
ostream& operator<<(ostream& os, const Vector2& vec) {
    os << "[" << vec.elements[0] << ", " << vec.elements[1] << ", " << vec.elements[2] << ", " << vec.elements[3] << "]";
    return os;
}
