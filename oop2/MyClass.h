#ifndef MYCLASS_H
#define MYCLASS_H

#include <iostream>
#include <fstream>
using namespace std;

class MyClass {
private:
    int a;
    double b;

public:
    // Конструкторы
    MyClass(); // Конструктор по умолчанию
    MyClass(int A, double B); // Конструктор с параметрами
    MyClass(const MyClass& other); // Конструктор копирования

    ~MyClass(); // Деструктор (сохранение данных в файл)

    void print() const; // Перегруженный метод печати
    void print(ostream& os) const; // Печать в указанный поток

    // Геттеры
    int getA() const;
    double getB() const;

    // Сеттеры (с проверкой корректности ввода)
    void setA(int A);
    void setB(double B);
};

#endif // MYCLASS_H
