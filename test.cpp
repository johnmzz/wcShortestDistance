#include <bits/stdc++.h>
using namespace std;

template<typename T>
bool compare(const T &a,const T &b){
    return a.second<b.second;
}

uint32_t findw(const vector<pair<uint16_t, uint8_t>>& vec, uint8_t r) {
    uint32_t st = 0;
    uint32_t end = vec.size() - 1;

    while (st <= end) {
        uint32_t mid = (st + end) / 2;
        if (vec[mid].second == r) {
            return mid;
        }
        else if (vec[mid].second < r) {
            st = mid + 1;
        }
        else {
            end = mid - 1;
        }
    }
    return st;
}

int main()
{
    // initializing the vector of pairs
    vector<pair<uint16_t, uint8_t>> vect;

    // insertion of pairs (key, value) in vector vect
    vect.push_back(make_pair(1, 20));
    vect.push_back(make_pair(3, 42));
    vect.push_back(make_pair(5, 66));

    // printing the sorted vector
    cout << "KEY" << '\t' << "ELEMENT" << endl;
    for (pair<uint16_t, uint8_t>& x : vect)
        cout << x.first << '\t' << x.second << endl;


    auto low = std::lower_bound(vect.begin(), vect.end(), make_pair(0,70),compare<pair<uint16_t,uint8_t>>);

    cout << "lower_bound = " << low - vect.begin() << endl;

    auto result = findw(vect, 70);

    cout << "find_w = " << result << endl;

    return 0;
}