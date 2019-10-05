#include <random>

using namespace std;


// base random variable generator
class BRV_generator{
private:
    random_device rd;
    mt19937 generator;
    uniform_real_distribution<double> uniform_distr;

public:
    BRV_generator(){
        generator = mt19937(rd());
        uniform_distr = uniform_real_distribution<double>(0.0, 1.0);
    }

    double roll_brv(){
        return uniform_distr(generator);
    }
};
