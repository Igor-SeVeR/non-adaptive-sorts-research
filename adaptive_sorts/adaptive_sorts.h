#ifndef COURSEWORK_ADAPTIVE_SORTS_H
#define COURSEWORK_ADAPTIVE_SORTS_H

#include <vector>
#include <ctime>
#include <random>
#include <iostream>

// Method to swap two numbers
template<typename T>
static void swap(T &a, T &b) {
    T c = a;
    a = b;
    b = c;
}

class heap_sort {
private:

    template<typename T>
    static void build_heap(std::vector<T> &arr) {
        size_t size = arr.size();
        for (int i = ((int)size) - 1; i >= 0; --i) {
            sift_down(arr, i, size);
        }
    }

    template<typename T>
    static void sift_down(std::vector<T> &arr, size_t pos, size_t cur_size) {
        size_t cur_elem = pos;
        while (cur_elem * 2 + 1 < cur_size) {
            size_t left_child_index = cur_elem * 2 + 1;
            size_t right_child_index = cur_elem * 2 + 2;
            size_t max_child_index = cur_elem;
            if (arr[left_child_index] > arr[max_child_index])
                max_child_index = left_child_index;
            if (right_child_index < cur_size && arr[right_child_index] > arr[max_child_index])
                max_child_index = right_child_index;
            if (arr[cur_elem] < arr[max_child_index]) {
                swap(arr[cur_elem], arr[max_child_index]);
                cur_elem = max_child_index;
            } else return;
        }
    }

public:

    template<typename T>
    static void sort(std::vector<T> &arr) {
        build_heap(arr);
        size_t size = arr.size();
        for (size_t i = 0; i < size - 1; ++i) {
            swap(arr[0], arr[size - i - 1]);
            sift_down(arr, 0, size - i - 1);
        }
    }

};

class merge_sort {
private:
    template<typename T>
    static void merge(std::vector<T> &arr, size_t l, size_t m, size_t r) {
        size_t index = 0;
        size_t buffer_size = r - l + 1;
        std::vector<T> buffer(buffer_size);
        size_t l_ind = l;
        size_t mid = m;
        while (index < buffer_size) {
            if (l_ind >= m)
                buffer[index++] = arr[mid++];
            else if (mid >= r)
                buffer[index++] = arr[l_ind++];
            else {
                if (arr[l_ind] <= arr[mid])
                    buffer[index++] = arr[l_ind++];
                else
                    buffer[index++] = arr[mid++];
            }
        }
        for (size_t i = l; i <= r; i++)
            arr[i] = buffer[i - l];
    }

    template<typename T>
    static void sort_merge(std::vector<T> &arr, size_t l, size_t r) {
        if (r <= l + 1)
            return;
        size_t m = (r + l) / 2;
        sort_merge(arr, l, m);
        sort_merge(arr, m, r);
        merge(arr, l, m, r);
    }

public:
    template<typename T>
    static void sort(std::vector<T> &arr) {
        sort_merge(arr, 0, arr.size());
    }

};


// Class, which implements adaptive heap sort in О(n log n) sorting + O(n) building
class qsort {
public:
    // Method to sort array using quick sort algorithm
    template<typename T>
    static void sort(std::vector<T> &arr, size_t l, size_t r) {
        srand(time(nullptr));
        if (r - l <= 1)
            return;
        T x = arr[l + rand() % (r - l)];
        int smaller_block = l;
        int equals_block = l;
        for (int i = l; i < r; ++i) {
            if (arr[i] < x) {
                swap(arr[i], arr[smaller_block]);
                if (smaller_block != equals_block)
                    swap(arr[i], arr[equals_block]);
                smaller_block++;
                equals_block++;
            }
            else if (arr[i] == x) {
                swap(arr[i], arr[equals_block]);
                equals_block++;
            }
        }
        sort(arr, l, smaller_block);
        sort(arr, equals_block, r);
    }

};

#endif //COURSEWORK_ADAPTIVE_SORTS_H
