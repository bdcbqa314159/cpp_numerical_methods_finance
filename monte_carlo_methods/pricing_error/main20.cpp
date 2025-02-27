#include <iostream>
#include "path_dep_option02.hpp"

using namespace std;

int main()
{
     double S0 = 100.0, r = 0.03, sigma = 0.2;
     BSModel Model(S0, r, sigma);

     double T = 1.0 / 12.0, K = 100.0;
     int m = 30;
     ArthmAsianCall Option(T, K, m);

     long N = 30000;
     cout << "Asian Call Price = "
          << Option.PriceByMC(Model, N) << endl
          << "Pricing Error = " << Option.PricingError << endl;

     return 0;
}
