#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

using namespace std;

enum ft_direction {
    forward,
    inverse
};

ostream& operator<<(ostream& out, const vector<complex<double>>& v) {
    for (auto& elem : v) {
        out << elem << " ";
    }
    out << endl;
    return out;
}

//discrete fourirer transform for a one-dimensional function
void dft1d(const vector<complex<double>>& in, vector<complex<double>>& out, ft_direction direction) {
    for (int idx_out = 0; idx_out < in.size(); ++idx_out) {
        complex<double> sum(0, 0);
        for (int idx_in = 0; idx_in < in.size(); ++idx_in) {
            double expnt = 2 * M_PI * idx_in * idx_out / in.size();
            expnt *= direction == ft_direction::inverse ? 1 : -1;
            sum += in[idx_in] * exp(complex<double>(0, expnt));
        }
        sum /= direction == ft_direction::inverse ? in.size() : 1;
        out.push_back(sum);
    }
}

int main(int argc, char* argv[]) {
    const int size = 256;

    vector<complex<double>> signal;
    double delta = 2 * M_PI / size;

    for (int i = 0; i < size; ++i) {
        double x = i * delta;
        signal.emplace_back(sin(x), 0); //add sample
    }

    cout << "original signal: " << endl;
    cout << signal;

    vector<complex<double>> transformed;
    dft1d(signal, transformed, ft_direction::forward);

    cout << "transformed signal: " << endl;
    cout << transformed;

    dft1d(transformed, signal, ft_direction::inverse);

    cout << "restored signal: " << endl;
    cout << signal;

    return 0;
}
