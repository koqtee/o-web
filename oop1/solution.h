#ifndef SOLUTION_H
#define SOLUTION_H

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

class Solution {
private:
    double a, b, c;

public:
    Solution(double A = 0, double B = 0, double C = 0); // Конструктор
    void print(); // Метод печати данных
    void solve(); // Метод вычисления корней
    string rootType(); // Метод определения типа корней
    void readFromFile(const string& filename); // Чтение коэффициентов из файла
};

#endif // SOLUTION_H
