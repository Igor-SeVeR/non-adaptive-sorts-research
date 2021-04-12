#include <iostream>
#include <vector>
#include <adaptive_sorts.h>
#include <nonadaptive_sorts.h>
#include "generator.h"

using namespace std;

int main() {
    generator gen;
    vector<int> arr = gen.generate_int(8, false, "");
    batchers_sort::sort(arr);
    for (int i = 0; i < arr.size(); i++)
        cout << arr[i] << '\n';
}