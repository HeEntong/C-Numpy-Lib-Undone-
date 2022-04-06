#ifndef HEAD_VEC
#define HEAD_VEC
#include <vector>
#include <cstdio>
#endif


template <typename T>
int partition(std::vector<T> &, int, int);

template <typename T>
void fastSort(std::vector<T> &, int, int);

template <typename T>
int partition(std::vector<T> &vec, int lower, int upper){
    int leftBound = lower;
    T pivot = vec[upper - 1];
    for (int i = lower; i < upper; i++){
        if (vec[i] < pivot){
            std::swap(vec[i], vec[leftBound]);
            leftBound++;
        }
    }
    std::swap(vec[upper - 1], vec[leftBound]);
    return leftBound;
}

template <typename T>
void fastSort(std::vector<T> &vec, int lower, int upper){
    if (upper - lower == 1 || upper == lower){
        return;
    }
    else{
        auto boundary = partition(vec, lower, upper);
        fastSort(vec, lower, boundary);
        fastSort(vec, boundary + 1, upper);
    }
} 