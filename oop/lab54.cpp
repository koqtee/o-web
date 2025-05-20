#include <iostream>
#include <string>
using namespace std;

int main() {
    string text;
    cout << "Введите сообщение: ";
    getline(cin, text);

    int maxCount = 0, currentCount = 0;

    for (char c : text) {
        if (isdigit(c)) {
            currentCount++;
            maxCount = max(maxCount, currentCount);
        } else {
            currentCount = 0;
        }
    }

    cout << "Наибольшее количество цифр, идущих подряд: " << maxCount << endl;
    return 0;
}
