#include <iostream>
#include <vector>
#include <random>

using namespace std;

int main() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist_1(1, 6);
    std::uniform_int_distribution<> dist_2(10,20);

    for (int n=0; n<10; n++) {
        std::cout << dist_1(gen) << ',';
        std::cout << dist_2(gen) << ' ';
    }
    std::cout << '\n';
}
