#include <iostream>
using namespace std;

int main() {
    int mas[15]; // Исходный массив
    int newMas[15]; // Новый массив
    int index = 0;

    // Ввод исходного массива
    cout << "Введите 15 элементов массива: ";
    for (int i = 0; i < 15; i++) {
        cin >> mas[i];
    }

    // Сначала добавляем элементы с чётных индексов
    for (int i = 0; i < 15; i += 2) {
        newMas[index++] = mas[i];
    }

    // Затем добавляем элементы с нечётных индексов
    for (int i = 1; i < 15; i += 2) {
        newMas[index++] = mas[i];
    }

    // Вывод нового массива
    cout << "Новый массив: ";
    for (int i = 0; i < 15; i++) {
        cout << newMas[i] << " ";
    }

    return 0;
}
