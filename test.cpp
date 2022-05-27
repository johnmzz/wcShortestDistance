#include <bits/stdc++.h>
using namespace std;

template<typename T>
bool compare(const T &a,const T &b){
    return a.second<b.second;
}

int main()
{
    // initializing the vector of pairs
    vector<pair<int, int> > vect;

    // insertion of pairs (key, value) in vector vect
    vect.push_back(make_pair(1, 20));
    vect.push_back(make_pair(3, 42));
    vect.push_back(make_pair(5, 66));

    // printing the sorted vector
    cout << "KEY" << '\t' << "ELEMENT" << endl;
    for (pair<int, int>& x : vect)
        cout << x.first << '\t' << x.second << endl;


    auto low = std::lower_bound(vect.begin(), vect.end(), make_pair(0,70),compare<pair<int,int>>);

    cout << "lower_bound = " << low - vect.begin() << endl;

    return 0;
}