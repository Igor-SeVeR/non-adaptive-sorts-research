//
// Created by igor_ on 13.04.2021.
//

#ifndef COURSEWORK_GENERATOR_H
#define COURSEWORK_GENERATOR_H

#include <random>
#include <time.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <climits>

class generator {
public:

    generator() : rnd(time(nullptr)) {
    }

    std::vector<int> generate_int(int number, bool need_to_write_to_file, const std::string &filename) {
        std::vector<int> arr(number);
        for (size_t i = 0; i < number; ++i) {
            arr[i] = rnd() % INT_MAX;
            if (rnd() % 2 == 1)
                arr[i] = -arr[i];
        }
        if (need_to_write_to_file) {
            std::freopen(filename.c_str(), "w", stdout);
            std::cout << number << '\n';
            for (size_t i = 0; i < number; ++i)
                std::cout << arr[i] << ' ';
            std::fclose(stdout);
        }
        return arr;
    }

private:
    std::mt19937 rnd;
};


#endif //COURSEWORK_GENERATOR_H
