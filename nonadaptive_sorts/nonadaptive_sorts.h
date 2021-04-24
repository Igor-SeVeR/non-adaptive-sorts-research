//
// Created by igor_ on 12.04.2021.
//

#ifndef COURSEWORK_NONADAPTIVE_SORTS_H
#define COURSEWORK_NONADAPTIVE_SORTS_H

#include <vector>
#include <iostream>
#include <omp.h>
#include <ctime>
#include <random>
#include <algorithm>

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

    template<typename T>
    static void merge_sort(std::vector<T> &arr, size_t l, size_t r, size_t real_size) {
        if (r <= l)
            return;
        size_t m = (r + l) / 2;
        merge_sort(arr, l, m, real_size);
        merge_sort(arr, m + 1, r, real_size);
        merge(arr, l, m, r, real_size);
    }

public:

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
    static void merge_if_50000(std::vector<T> &arr, size_t l, size_t m, size_t r, size_t real_size, int number_of_threads) {
        omp_set_dynamic(0);
        omp_set_num_threads(number_of_threads);
        size_t buffer_size = r - l + 1;
        size_t i, k = buffer_size / 2, j;
        for (k; k > 100; k /= 2) {
#pragma omp for private(i, j)
                for (j = k % (buffer_size / 2); j < buffer_size - k; j += k + k) {
                    for (i = 0; i < k; i++)
                        if (l + j + i < real_size && l + j + i + k < real_size)
                            operations::compexch(arr[l + j + i], arr[l + j + i + k]);
                }
            }
        for (k; k > 0; k /= 2) {
            for (j = k % (buffer_size / 2); j < buffer_size - k; j += k + k) {
                for (i = 0; i < k; i++) {
                    if (l + j + i < real_size && l + j + i + k < real_size)
                        operations::compexch(arr[l + j + i], arr[l + j + i + k]);
                }
            }
        }
    }

    template<typename T>
    static void merge(std::vector<T> &arr, size_t l, size_t m, size_t r, size_t real_size, int number_of_threads) {
        size_t buffer_size = r - l + 1;
        size_t i, k, j;
        for (k = buffer_size / 2; k > 0; k /= 2) {
            for (j = k % (buffer_size / 2); j < buffer_size - k; j += k + k) {
                for (i = 0; i < k; i++) {
                    if (l + j + i < real_size && l + j + i + k < real_size)
                        operations::compexch(arr[l + j + i], arr[l + j + i + k]);
                }
            }
        }
    }

    template<typename T>
    static void merge_sort(std::vector<T> &arr, size_t l, size_t r, size_t real_size, int number_of_threads) {
        if (r <= l)
            return;
        size_t m = (r + l) / 2;
        merge_sort(arr, l, m, real_size, number_of_threads);
        merge_sort(arr, m + 1, r, real_size, number_of_threads);
        //if (r - l + 1 >= 10000)
            merge_if_50000(arr, l, m, r, real_size, number_of_threads);
        //else
        //    merge(arr, l, m, r, real_size, number_of_threads);

    }

public:

    template<typename T>
    static void sort(std::vector<T> &arr, int number_of_threads = 8) {

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
        merge_sort(arr, 0, closest_right_st - 1, real_size, number_of_threads);
    }

};

class zig_zag_sort {

private:

    template<class RandomType>
    static inline std::vector<int> gen_perm(int n, RandomType& rnd) {
        std::vector<int> p(n);
        std::iota(p.begin(), p.end(), 0);
        std::shuffle(p.begin(), p.end(), rnd);
        return p;
    }

    template<class RandomType>
    static inline void epsHalver(int l1, int r1, int l2, int r2, int rev_eps, RandomType& rnd, std::vector<int> &arr, int real_size) {
        int n = r1 - l1;
        int r = std::max((int)std::ceil(2 * rev_eps * log(rev_eps)), 1);
        for (int iter = 0; iter < r; ++iter) {
            std::vector<int> p = gen_perm(n, rnd);
            for (int i = 0; i < n; ++i) {
                if (i + l1 < real_size && p[i] + l2 < real_size)
                    operations::compexch(arr[i + l1], arr[p[i] + l2]);
            }
        }
    }

    static inline void bubbleSort(int l1, int r1, int l2, int r2, std::vector<int> &arr, int real_size) {
        auto get_idx = [&](int i) {
            if (i < r1 - l1) {
                return l1 + i;
            }
            i -= r1 - l1;
            return l2 + i;
        };
        int n = r1 - l1;
        n <<= 1;
        for (int i = 0; i < n - 1; ++i) {
            for (int j = 0; j < n - 1 - i; ++j) {
                if (get_idx(j + 1) < real_size)
                    operations::compexch(arr[get_idx(j)], arr[get_idx(j + 1)]);
            }
        }
    }

    template<class RandomType>
    static void Attenuate(int l1, int r1, int l2, int r2, int rev_eps, RandomType& rnd, std::vector<int> &arr, int real_size) {
        int n = r1 - l1;
        if (n <= 8) {
            bubbleSort(l1, r1, l2, r2, arr, real_size);
            return;
        }
        int m1 = (l1 + r1) >> 1, m2 = (l2 + r2) >> 1;

        epsHalver(l1, m1, m1, r1, rev_eps, rnd, arr, real_size);
        epsHalver(l2, m2, m2, r2, rev_eps, rnd, arr, real_size);
        epsHalver(m1, r1, l2, m2, rev_eps, rnd, arr, real_size);
        Attenuate(m1, r1, l2, m2, rev_eps, rnd, arr, real_size);

        l1 = m1; r2 = m2;
        m1 = (l1 + r1) >> 1, m2 = (l2 + r2) >> 1;

        epsHalver(l1, m1, m1, r1, rev_eps, rnd, arr, real_size);
        epsHalver(l2, m2, m2, r2, rev_eps, rnd, arr, real_size);
        epsHalver(m1, r1, l2, m2, rev_eps, rnd, arr, real_size);
        Attenuate(m1, r1, l2, m2, rev_eps, rnd, arr, real_size);
    }

    template<class RandomType>
    static void Reduce(int l1, int r1, int l2, int r2, int rev_eps, RandomType& rnd, std::vector<int> &arr, int real_size) {
        int n = r1 - l1;
        if (n <= 8) {
            bubbleSort(l1, r1, l2, r2, arr, real_size);
            return;
        }
        epsHalver(l1, r1, l2, r2, rev_eps, rnd, arr, real_size);
        Attenuate(l1, r1, l2, r2, rev_eps, rnd, arr, real_size);
    }

public:

    static void sort(std::vector<int> &arr) {

        size_t closest_right_st = 0;
        size_t closest_right_st_val = 1;
        size_t size = arr.size();
        size_t real_size = size;
        int eps = 15;
        std::mt19937 rnd(time(nullptr));
        --size;
        while (size > 0) {
            size = size >> 1;
            closest_right_st = closest_right_st += 1;
            closest_right_st_val <<= 1;
        }
        size = real_size;
        for (int j = 1; j <= closest_right_st; ++j) {
            size_t part_size = closest_right_st_val / (1 << (j));
            for (int i = 1; i <= 1 << (j - 1); ++i) {
                Reduce((i - 1) * (2 * part_size), (i - 1) * (2 * part_size) + part_size, i * 2 * part_size - part_size, i * 2 * part_size, eps, rnd, arr, size);
            }
            for (int i = 1; i <= (1 << (j)) - 1; i++) {
                for (int k = 0; k < part_size; k++)
                    if (i * part_size + k < size)
                        operations::compexch(arr[(i - 1) * part_size + k], arr[i * part_size + k]);
                Reduce((i - 1) * part_size, i * part_size, i * part_size, (i + 1) * part_size, eps, rnd, arr, size);
            }
            for (int i = (1 << j); i >= 2; i--) {
                for (int k = 0; k < part_size; k++)
                    if ((i - 1) * part_size + k < size)
                        operations::compexch(arr[(i - 1) * part_size + k], arr[(i - 2) * part_size + k]);
                Reduce((i - 2) * part_size, (i - 1) * part_size, (i - 1) * part_size, i * part_size, eps, rnd, arr, size);
            }
        }
    }
};

#endif //COURSEWORK_NONADAPTIVE_SORTS_H
