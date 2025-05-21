#include "Solution.h"

int main() {
    Solution eq;
    eq.readFromFile("input.txt"); // Ввод коэффициентов из файла
    eq.print(); // Вывод коэффициентов
    eq.solve(); // Решение уравнения
    cout << "Тип корней: " << eq.rootType() << endl;

    // Создание массива объектов
    Solution equations[2] = {Solution(1, -3, 2), Solution(1, 2, 1)};
    
    cout << "\nМассив уравнений:" << endl;
    for (int i = 0; i < 2; i++) {
        equations[i].print();
        equations[i].solve();
    }

    return 0;
}
