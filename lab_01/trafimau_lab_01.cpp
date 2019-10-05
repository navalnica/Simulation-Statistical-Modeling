/*
 * This program creates two pseudo-random generators
 * with uniform distribution on [0; 1].
 * The first one uses multiplicative congruential method.
 * The second one combines the first generator with the
 * c++ builtin uniform random generator to use MacLaren-Marsaglia method.
 * Results are tested with Chi-squared test and Kolmogorov test.
 * */


#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <climits>
#include <numeric>
#include <random>


using namespace std;

using ull = unsigned long long;

// parameters
ull const A_STAR_0 = 24'389;
ull const BETA = A_STAR_0;
ull const M = 1ULL << 31ULL;
int const K = 32;
double const chi_square_threshold = 16.92;
double const kolmogorov_threshold = 1.36;

vector<double> mult_congr(int n = 1000) {
    // multiplicative congruent(ial generator
    vector<ull> a_star(n, -111);
    a_star[0] = A_STAR_0;
    for (int i = 1; i < n; ++i) {
        a_star[i] = (BETA * a_star[i - 1]) % M;
    }
    vector<double> a(n, -111);
    for (int i = 0; i < n; ++i) {
        a[i] = static_cast<double>(a_star[i]) / M;
    }
    return a;
}

vector<double> builtin_uniform(int n = 1000) {
    vector<double> res(n);
    random_device rd; // obtain a seed
    mt19937 engine(rd()); // create engine and seed it with rd
    uniform_real_distribution<double> distribution(
            0.0,
            nextafter(1.0, numeric_limits<double>::max())
    );
    for (int i = 0; i < n; ++i) {
        res[i] = distribution(engine);
    }
    return res;
}

vector<double> maclaren_marsaglia(
        vector<double> &first,
        vector<double> &second
) {
    int size = min(first.size(), second.size()) - K;
    vector<double> res(size);
    vector<double> v(first.begin(), first.begin() + K);
    for (int i = 0; i < size; ++i) {
        int s = int(second[i] * K);
        res[i] = v[s];
        v[s] = first[i + K];
    }
    return res;
}

void simple_tests(vector<double> &sample) {
    cout << "size: " << sample.size() << endl;
    cout << "some elements:" << endl;
    for (int i = 0; i < 5; ++i) {
        cout << sample[i * 10] << " ";
    }
    cout << endl;

    double mean = accumulate(
            sample.begin(),
            sample.end(),
            0.0
    ) / sample.size();
    cout << "mean: " << mean << endl;

    double std = 0;
    for (double d : sample) {
        std += (d - mean) * (d - mean);
    }
    std /= (sample.size() - 1);
    cout << "std: " << std << endl;
}

vector<int> calc_bins(
        vector<double> &v, int bins_cnt = 10,
        double start = 0.0, double end = 1.0
) {
    // calculate histogram

    // copy vector to sort inplace
    vector<double> copy(v.cbegin(), v.cend());
    sort(copy.begin(), copy.end());

    vector<int> res(bins_cnt);
    const double stride = (end - start) / bins_cnt;
    int prev = 0;
    for (int i = 0; i < bins_cnt; ++i) {
        double threshold = stride * (i + 1);
        auto it = upper_bound(copy.cbegin(), copy.cend(), threshold); // min greater element's iterator
        int cnt = distance(copy.cbegin(), it);
        res[i] = cnt - prev;
        prev = cnt;
    }
    return res;
}

double chi_square(vector<double> &v) {
    // check distribution against uniform on [0; 1]

    int n = v.size();
    vector<int> bins = calc_bins(v);
    double res = 0;
    double prob = n / 10.0; // P(a <= x <= a + 1/10)
    for (auto cnt : bins) {
        res += (cnt - prob) * (cnt - prob) / prob;
    }
    return res;
}

double kolmogorov(vector<double> &v) {
    // check distribution against uniform on [0; 1]

    // copy vector to sort inplace
    vector<double> copy(v.cbegin(), v.cend());
    sort(copy.begin(), copy.end());

    double res = 0;
    for (int i = 0; i < copy.size(); ++i) {
        double const f = min(max(copy[i], 0.0), 1.0); // CDF value of the reference distribution
        double const f_hat = static_cast<double>(i + 1) / copy.size(); // empirical CDF value
        double cur = abs(f_hat - f);
        res = max(res, cur);
    }
    return res;
}

void main_tests(vector<double> &v, string const &method_name) {
    cout << endl << "chi-square test for " << method_name << ":" << endl;
    double chi_square_val = chi_square(v);
    cout << "value: " << chi_square_val << endl;
    cout << "test passed: " << boolalpha << (chi_square_val < chi_square_threshold) << endl;

    cout << endl << "kolmogorov test for " << method_name << ":" << endl;
    double kolmogorov_value = kolmogorov(v);
    cout << "value: " << kolmogorov_value << endl;
    cout << "test passed: " << boolalpha <<
         (sqrt(v.size()) * kolmogorov_value < kolmogorov_threshold) << endl;
}

int main() {
    string const separator = "---------------------------";

    cout << endl << separator << endl << endl;
    cout << "parameters:" << endl;
    cout << "A_STAR_0: " << A_STAR_0 << endl;
    cout << "M: " << M << endl;
    cout << "K: " << K << endl;

    cout << endl << separator << endl << endl;

    cout << "multiplicative congruential generator:" << endl;
    vector<double> mult_congr_res = mult_congr();
    simple_tests(mult_congr_res);
    main_tests(mult_congr_res, "multiplicative congruential");

    cout << endl << separator << endl << endl;

    cout << endl << "builtin uniform generator:" << endl;
    vector<double> builtin_uniform_res = builtin_uniform();
    simple_tests(builtin_uniform_res);
    main_tests(builtin_uniform_res, "builtin uniform");

    cout << endl << separator << endl << endl;

    cout << endl << "maclaren marsaglia:" << endl;
    vector<double> maclaren_marsaglia_res = maclaren_marsaglia(mult_congr_res, builtin_uniform_res);
    simple_tests(maclaren_marsaglia_res);
    main_tests(maclaren_marsaglia_res, "maclaren marsaglia");

    return 0;
}