#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
using namespace std;

int main() {
    const int rows = 7;
    const int cols = 7;
    int P;
    int W[rows][cols], Z[rows][cols];

    srand(time(0)); // Инициализируем генератор случайных чисел

    // Заполняем матрицу W случайными числами от 0 до 9
    cout << "Матрица W:" << endl;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            W[i][j] = rand() % 10;
            cout << W[i][j] << " ";
        }
        cout << endl;
    }

    // Заполняем матрицу Z случайными числами от 0 до 9
    cout << "\nМатрица Z:" << endl;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            Z[i][j] = rand() % 10;
            cout << Z[i][j] << " ";
        }
        cout << endl;
    }

    cout << "\nВведите число P: ";
    cin >> P;

    vector<int> T, S;

    // Формируем массив T из элементов матрицы W, больших P
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (W[i][j] > P) {
                T.push_back(W[i][j]);
            }
        }
    }

    // Формируем массив S из элементов матрицы Z, больших P
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (Z[i][j] > P) {
                S.push_back(Z[i][j]);
            }
        }
    }

    // Вывод массивов T и S
    cout << "\nМассив T (элементы W > P): ";
    for (int num : T) {
        cout << num << " ";
    }
    cout << endl;

    cout << "Массив S (элементы Z > P): ";
    for (int num : S) {
        cout << num << " ";
    }
    cout << endl;

    return 0;
}
