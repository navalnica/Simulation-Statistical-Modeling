#include <cmath>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <iomanip>

#include "Moments.h"
#include "BRV_generator.h"


using namespace std;


vector<int> get_bernoulli(double p, int size){
    vector<int> res;
    BRV_generator brv_generator;
    for(int i = 0; i < size; ++i){
        double brv = brv_generator.roll_brv();
        int roll = int(brv <= p);
        res.push_back(roll);
    }
    return res;
}


vector<int> get_binomial(int m, double p, int size){
    vector<int> res;
    BRV_generator brv_generator;
    for(int i = 0; i < size; ++i){
        int sum = 0;
        for(int j = 0; j < m; ++j){
            double brv = brv_generator.roll_brv();
            sum += int(brv <= p);
        }
        res.push_back(sum);
    }
    return res;
}

vector<int> get_geometric(double p, int size){
    vector<int> res(size, 0);
    BRV_generator brv_generator;
    for(int i = 0; i < size; ++i){
        double brv = brv_generator.roll_brv();
        double tmp = log(brv) / log(1 - p);
        int roll = int(ceil(tmp));
        res[i] = roll;
    }
    return res;
}


vector<int> get_poisson(int lambda, int size){
    vector<int> res(size, 0);
    BRV_generator brv_generator;
    for(int i = 0; i < size; ++i){
        int k = 1;
        double bsv = brv_generator.roll_brv();
        double product = bsv;
        while (product >= exp(-lambda)){
            bsv = brv_generator.roll_brv();
            product *= bsv;
            k++;
        }
        res[i] = k - 1;
    }
    return res;
}


void perform_tests(vector<int> const & vec){
    Moments moments(vec);
    cout << "mean: " << moments.mean() << endl;
    cout << "variance: " << moments.variance_unbiased() << endl;
    cout << "skewness: " << moments.skewness_unbiased() << endl;
    cout << "kurtosis: " << moments.kurtosis_unbiased() << endl;
}

void write_to_file(vector<int> const & vec, ofstream & out){
    out << fixed << setprecision(3);
    for(auto el : vec) { out << el << endl; }
}

int main(){
    string separator = "--------------------";

    cout << fixed << setprecision(3);
    int const n_rolls = 1000;

    
    // -- bernoulli --
    
    cout << separator << endl;
    cout << "bernoulli:" << endl;
    
    double const bernoulli_p = 0.7;
    auto bernoulli = get_bernoulli(bernoulli_p, n_rolls);
    perform_tests(bernoulli);

    ostringstream bernoulli_filename;
    bernoulli_filename << fixed << setprecision(2);
    bernoulli_filename << "bernoulli_" << bernoulli_p << ".txt";
    ofstream bernoulli_out(bernoulli_filename.str());
    write_to_file(bernoulli, bernoulli_out);
    


    // -- binomial --
    
    cout << endl << separator << endl;
    cout << "binomial" << endl;
    
    int const binomial_m = 5;
    double const binomial_p = 0.25;
    auto binomial = get_binomial(binomial_m, binomial_p, n_rolls);
    perform_tests(binomial);
    
    ostringstream binomial_filename;
    binomial_filename << fixed << setprecision(2);
    binomial_filename << "binomial_" << binomial_m << "_" << binomial_p << ".txt";
    ofstream binomial_out(binomial_filename.str());
    write_to_file(binomial, binomial_out);


    // -- geometric --

    cout << endl << separator << endl;
    cout << "geometric:" << endl;

    double geometric_p = 0.7;
    auto geometric = get_geometric(geometric_p, n_rolls);
    perform_tests(geometric);

    ostringstream geometric_filename;
    geometric_filename << fixed << setprecision(2);
    geometric_filename << "geometric_" << geometric_p << ".txt";
    ofstream geometric_out(geometric_filename.str());
    write_to_file(geometric, geometric_out);


    // -- poisson--

    cout << endl << separator << endl;
    cout << "poisson:" << endl;

    int poisson_lambda = 2;
    auto poisson = get_poisson(poisson_lambda, n_rolls);
    perform_tests(poisson);

    ostringstream poisson_filename;
    poisson_filename << fixed << setprecision(2);
    poisson_filename << "poisson_" << poisson_lambda << ".txt";
    ofstream poisson_out(poisson_filename.str());
    write_to_file(poisson, poisson_out);

    return 0;
}
