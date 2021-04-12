//
// Created by igor_ on 12.04.2021.
//

#ifndef COURSEWORK_NONADAPTIVE_SORTS_H
#define COURSEWORK_NONADAPTIVE_SORTS_H

#include <vector>

// Class with all operations needed
class operations {
public:

    // Method to swap two numbers
    template<typename T>
    static void swap(T &a, T &b) {
        T c = a;
        a = b;
        b = c;
    }

    // Method to compare and excange two T elements
    template<typename T>
    static void compexch(T &a, T &b) {
        if (a > b)
            swap(a, b);
    }
};

class batchers_sort {

private:
    template<typename T>
    static void merge(std::vector<T> &arr, size_t l, size_t m, size_t r) {
        size_t buffer_size = r - l + 1;
        for (size_t k = buffer_size / 2; k > 0; k /= 2)
            for (size_t j = k % (buffer_size / 2); j + k < buffer_size; j += k + k)
                for (size_t i = 0; i < k; i++)
                    operations::compexch(arr[l + j + i], arr[l + j + i + k]);
    }

public:
    template<typename T>
    static void merge_sort(std::vector<T> &arr, size_t l, size_t r) {
        if (r <= l)
            return;
        size_t m = (r + l) / 2;
        merge_sort(arr, l, m);
        merge_sort(arr, m + 1, r);
        merge(arr, l, m, r);
    }

    template<typename T>
    static void sort(std::vector<T> &arr) {
        merge_sort(arr, 0, arr.size() - 1);
    }

};

class bubble_sort {
public:

    //Method to sort algorithm using bubble sort in O(n ^ 2)
    template<typename T>
    static void sort(std::vector<T> &arr) {\
    size_t size = arr.size();
        for (size_t i = 0; i < size; ++i) {
            for (size_t j = 0; j < size - i - 1; ++j) {
                operations::compexch(arr[j], arr[j + 1]);
            }
        }
    }
};


class nonadaptive_heap_sort {
private:

    template<typename T>
    static void sift_up(std::vector<T> &arr, size_t pos) {
        while (pos > 0) {
            size_t parents_index = (pos - 1) / 2;
            operations::compexch(arr[pos], arr[parents_index]);
            pos = parents_index;
        }
    }

public:
    template<typename T>
    static void sort(std::vector<T> &arr) {
        size_t size = arr.size();
        for (size_t i = 0; i < size - 1; ++i) {
            for (size_t j = size - i - 1; j > 0; --j)
                sift_up(arr, j);
            operations::swap(arr[0], arr[size - i - 1]);
        }
    }
};


#endif //COURSEWORK_NONADAPTIVE_SORTS_H