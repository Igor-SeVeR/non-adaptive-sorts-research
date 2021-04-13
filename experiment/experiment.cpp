#include <iostream>
#include <vector>
#include <adaptive_sorts.h>
#include <nonadaptive_sorts.h>
#include <chrono>
#include "generator.h"

using namespace std;

int main() {
    generator gen;
    vector<int> arr = gen.generate_int(1000000, false, "");
    vector<int> kek(arr.begin(), arr.end());
    unsigned int start = clock();
    qsort::sort(arr, 0, arr.size());
    unsigned int end = clock();
    for (int i = 0; i < arr.size() - 1; i++)
        if (arr[i] > arr[i + 1]) {
            std::cout << "Incorrect_sort";
            return 0;
        }
    std::cout << "Correct_sort" << ' ' << end - start << '\n';
    start = clock();
    batchers_parallel_sort::sort(kek, 8);
    end = clock();
    for (int i = 0; i < kek.size() - 1; i++)
        if (kek[i] > kek[i + 1]) {
            std::cout << "Incorrect_sort";
            return 0;
        }
    std::cout << "Correct_sort" << ' ' << end - start << '\n';
}