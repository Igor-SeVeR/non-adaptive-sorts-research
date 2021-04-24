#include <iostream>
#include <vector>
#include <adaptive_sorts.h>
#include <nonadaptive_sorts.h>
#include <chrono>
#include "generator.h"

using namespace std;

int main() {
    generator gen;
    vector<int> arr = gen.generate_int(1000, false, "");
    int n = 100;
    int sum = 0;
    for (int i = 0; i < n; i++) {
        vector<int> kek(arr.begin(), arr.end());
        unsigned int start = clock();
        qsort::sort(kek, 0, arr.size());
        unsigned int end = clock();
        sum += end - start;
    }
    cout << "qsort = " << sum / n << '\n';

    //testing merge
    sum = 0;
    for (int i = 0; i < n; i++) {
        vector<int> kek(arr.begin(), arr.end());
        unsigned int start = clock();
        merge_sort::sort(kek);
        unsigned int end = clock();
        sum += end - start;
    }
    cout << "merge = " << sum / n << '\n';

    //testing heap
    sum = 0;
    for (int i = 0; i < n; i++) {
        vector<int> kek(arr.begin(), arr.end());
        unsigned int start = clock();
        heap_sort::sort(kek);
        unsigned int end = clock();
        sum += end - start;
    }
    cout << "heap = " << sum / n << '\n';

    //testing bubble
    sum = 0;
    for (int i = 0; i < n; i++) {
        vector<int> kek(arr.begin(), arr.end());
        unsigned int start = clock();
        bubble_sort::sort(arr);
        unsigned int end = clock();
        sum += end - start;
    }
    cout << "bubble = " << sum / n << '\n';

    //testing non_par_batch
    sum = 0;
    for (int i = 0; i < n; i++) {
        vector<int> kek(arr.begin(), arr.end());
        unsigned int start = clock();
        batchers_sort::sort(kek);
        unsigned int end = clock();
        sum += end - start;
    }
    cout << "non_par_batch = " << sum / n << '\n';
    //testing par_batch
    //sum = 0;
    //for (int i = 0; i < n; i++) {
    //    vector<int> kek(arr.begin(), arr.end());
    //    unsigned int start = clock();
    //    batchers_parallel_sort::sort(kek, 16);
    //    unsigned int end = clock();
    //    sum += end - start;
    //}
    //cout << "par_batch = " << sum / n << '\n';

    //testing zig-zag
    sum = 0;
    for (int i = 0; i < n; i++) {
        vector<int> kek(arr.begin(), arr.end());
        unsigned int start = clock();
        zig_zag_sort::sort(arr);
        unsigned int end = clock();
        sum += end - start;
    }
    cout << "zig_zag = " << sum / n << '\n';
}