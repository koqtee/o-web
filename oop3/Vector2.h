#ifndef VECTOR2_H
#define VECTOR2_H

#include <iostream>
using namespace std;

class Vector2 {
private:
    int elements[4]; // Четыре числа в векторе

public:
    Vector2(); // Конструктор по умолчанию
    Vector2(int e1, int e2, int e3, int e4); // Конструктор с параметрами

    void print() const; // Метод печати

    // Перегрузка операторов
    Vector2 operator+(const Vector2& other) const; // Сложение двух векторов
    Vector2 operator*(int scalar) const; // Умножение на число
    bool operator==(const Vector2& other) const; // Оператор сравнения
    Vector2 operator-() const; // Унарный оператор минус (изменение знака всех элементов)

    // Перегрузка оператора вывода
    friend ostream& operator<<(ostream& os, const Vector2& vec);
};

#endif // VECTOR2_H
