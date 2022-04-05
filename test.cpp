#include <iostream>
#include <vector>

using namespace std;

int main() {
    vector<uint32_t> vec {10,20,30,30,40,50};

    for (auto v : vec) {
        cout << v << endl;
    }

    vector<uint32_t>::iterator low = lower_bound(vec.begin(), vec.end(), 30);

    cout << "lower_bound = " << (low - vec.begin()) << endl;
    cout << "vec.begin() = " << vec.begin() - vec.begin() << endl;
    cout << "vec.end() = " << vec.end() - vec.begin() << endl;
}
