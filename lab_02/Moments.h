#include <vector>
#include <cmath>

using namespace std;


class Moments{

private:
    vector<int> const & elements;

public:

    Moments(vector<int> const & elements): elements(elements) {}

    double mean(){
        double mean = 0;
        for(auto el : elements){ mean += el; }
        mean /= elements.size();
        return mean;
    }

    double central_moment(int order){
        double res = 0;
        double _mean = mean();
        for(auto el : elements) { res += pow(el - _mean , order); }
        return res / elements.size();
    }

    double central_moment_2_unbiased(){
        int m = elements.size();
        double cm2 = central_moment(2);
        double res = cm2 * m / (m - 1);
        return res;
    }

    double variance_unbiased() { return central_moment_2_unbiased(); }

    double central_moment_3_unbiased(){
        int m = elements.size();
        double cm3 = central_moment(3);
        double res = cm3 * m / (m - 1) * m / (m - 2);
        return res;
    }

    double skewness_unbiased(){
        double cm3u = central_moment_3_unbiased();
        double cm2u = central_moment_2_unbiased();
        double res = cm3u / pow(cm2u, 1.5);
        return res;
    }

    double kurtosis_unbiased() {
        size_t m = elements.size();
        double cm2 = central_moment(2);
        double cm4 = central_moment(4);
        double res = cm4 / cm2 / cm2 - 3 + 6.0 / (m + 1);
        res *= double(m - 1) / (m - 2);
        res *= double(m + 1) / (m - 3);
        return res;
    }
};
