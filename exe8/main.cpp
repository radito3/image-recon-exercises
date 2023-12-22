#include <iostream>
#include <cmath>
#include <cassert>
#include <gsl/gsl_wavelet.h>
#include <gsl/gsl_sort.h>

const int DATA_SIZE = 256;
const int NUM_COEF = 20;

void print(const double* signal, int size) {
    for (auto i = 0; i < size; i++) {
        std::cout << signal[i] << ' ';
    }
    std::cout << std::endl;
}

double* generate(int size) {
    assert(size <= DATA_SIZE);

    double delta = 2.0 * M_PI / size;
    double* result = new double[size];

    for (int i = 0; i < size; i++) {
        result[i] = sin(i * delta);
    }
    return result;
}

void dwave4(double* signal, int size) {
    gsl_wavelet* w = gsl_wavelet_alloc(gsl_wavelet_daubechies, 4);
    gsl_wavelet_workspace* work = gsl_wavelet_workspace_alloc(DATA_SIZE);

    gsl_wavelet_transform_forward(w, signal, 1, size, work);

    double* abscoeff = new double[size];
    for (int i = 0; i < size; i++) {
        abscoeff[i] = fabs(signal[i]);
    }

    size_t* indices = new size_t[size];
    gsl_sort_index(indices, abscoeff, 1, size);

    for (int i = 0; i + NUM_COEF < size; i++) {
        signal[indices[i]] = 0;
    }

    gsl_wavelet_transform_inverse(w, signal, 1, size, work);

    delete[] indices;
    delete[] abscoeff;
    gsl_wavelet_free(w);
    gsl_wavelet_workspace_free(work);
}

int main() {
    double* input = generate(DATA_SIZE);

    std::cout << "original signal: ";
    print(input, DATA_SIZE);

    dwave4(input, DATA_SIZE);

    std::cout << std::endl << "compressed and restored signal: ";
    print(input, DATA_SIZE);

    delete[] input;
    return 0;
}