//
// Created by igor_ on 12.04.2021.
//

#ifndef COURSEWORK_NONADAPTIVE_SORTS_H
#define COURSEWORK_NONADAPTIVE_SORTS_H

#include <vector>
#include <iostream>
#include <omp.h>

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
    static void merge(std::vector<T> &arr, size_t l, size_t m, size_t r, size_t real_size) {
        size_t buffer_size = r - l + 1;
        for (size_t k = buffer_size / 2; k > 0; k /= 2)
            for (size_t j = k % (buffer_size / 2); j + k < buffer_size; j += k + k)
                for (size_t i = 0; i < k; i++)
                    if (l + j + i < real_size && l + j + i + k < real_size)
                        operations::compexch(arr[l + j + i], arr[l + j + i + k]);
    }

public:
    template<typename T>
    static void merge_sort(std::vector<T> &arr, size_t l, size_t r, size_t real_size) {
        if (r <= l)
            return;
        size_t m = (r + l) / 2;
        merge_sort(arr, l, m, real_size);
        merge_sort(arr, m + 1, r, real_size);
        merge(arr, l, m, r, real_size);
    }

    template<typename T>
    static void sort(std::vector<T> &arr) {
        size_t closest_right_st = 1;
        size_t size = arr.size();
        size_t real_size = size;
        if (size == 0)
            return;
        --size;
        while (size > 0) {
            size = size >> 1;
            closest_right_st = closest_right_st << 1;
        }
        merge_sort(arr, 0, closest_right_st - 1, real_size);
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

class batchers_parallel_sort {

private:
    template<typename T>
    static void merge(std::vector<T> &arr, size_t l, size_t m, size_t r, size_t real_size) {
        size_t buffer_size = r - l + 1;
        omp_set_dynamic(0);
        omp_set_num_threads(8);
        for (size_t k = buffer_size / 2; k > 0; k /= 2)
            for (size_t j = k % (buffer_size / 2); j + k < buffer_size; j += k + k) {
                size_t i;
                #pragma omp parallel for shared(arr, l, j, k, real_size) private(i)
                for (i = 0; i < k; i++)
                    if (l + j + i < real_size && l + j + i + k < real_size)
                        operations::compexch(arr[l + j + i], arr[l + j + i + k]);
            }
    }

public:
    template<typename T>
    static void merge_sort(std::vector<T> &arr, size_t l, size_t r, size_t real_size) {
        if (r <= l)
            return;
        size_t m = (r + l) / 2;
        merge_sort(arr, l, m, real_size);
        merge_sort(arr, m + 1, r, real_size);
        merge(arr, l, m, r, real_size);
    }

    template<typename T>
    static void sort(std::vector<T> &arr) {
        size_t closest_right_st = 1;
        size_t size = arr.size();
        size_t real_size = size;
        if (size == 0)
            return;
        --size;
        while (size > 0) {
            size = size >> 1;
            closest_right_st = closest_right_st << 1;
        }
        merge_sort(arr, 0, closest_right_st - 1, real_size);
    }

};


#endif //COURSEWORK_NONADAPTIVE_SORTS_H
