#include <iostream>
#include <sstream>
#include <map>
using namespace std;

int main() {
    string text;
    int p;

    cout << "Введите сообщение: ";
    getline(cin, text);
    cout << "Введите минимальное количество повторов (p): ";
    cin >> p;

    map<string, int> wordCount;
    stringstream ss(text);
    string word;

    while (ss >> word) {
        wordCount[word]++;
    }

    cout << "Слова, встречающиеся более " << p << " раз:" << endl;
    for (const auto& pair : wordCount) {
        if (pair.second > p) {
            cout << pair.first << " (" << pair.second << " раз(а))" << endl;
        }
    }
    
    // cout << pair.first << endl;
    
    return 0;
}

