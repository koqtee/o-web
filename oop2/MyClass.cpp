#include "MyClass.h"

// Конструктор по умолчанию
MyClass::MyClass() : a(0), b(0.0) {}

// Конструктор с параметрами
MyClass::MyClass(int A, double B) {
    setA(A);
    setB(B);
}

// Конструктор копирования
MyClass::MyClass(const MyClass& other) : a(other.a), b(other.b) {}

// Деструктор (сохранение данных в файл)
MyClass::~MyClass() {
    ofstream file("output.txt", ios::app);
    if (file.is_open()) {
        file << "Удалён объект с данными: a = " << a << ", b = " << b << endl;
        file.close();
    }
}

// Метод печати
void MyClass::print() const {
    cout << "a = " << a << ", b = " << b << endl;
}

// Перегруженный метод печати в поток
void MyClass::print(ostream& os) const {
    os << "a = " << a << ", b = " << b << endl;
}

// Геттеры
int MyClass::getA() const { return a; }
double MyClass::getB() const { return b; }

// Сеттеры (с проверкой корректности)
void MyClass::setA(int A) {
    if (A >= 0) {
        a = A;
    } else {
        cout << "Ошибка: значение a не может быть отрицательным!" << endl;
    }
}

void MyClass::setB(double B) {
    if (B >= 0.0) {
        b = B;
    } else {
        cout << "Ошибка: значение b не может быть отрицательным!" << endl;
    }
}
