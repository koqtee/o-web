#include <iostream>
#include <string>
using namespace std;

int main() {
    string str;
    cout << "Введите строку: ";
    getline(cin, str);

    int count = 0;
    for (char c : str) {
        if (c >= 'a' && c <= 'z') {
            count++;
        }
    }

    cout << "Количество прописных латинских букв: " << count << endl;
    return 0;
}