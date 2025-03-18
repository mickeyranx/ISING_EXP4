#include "test.h"
#include <random>



uniform_real_distribution<double> test::distr = uniform_real_distribution<double>(0.0, 1.0);
default_random_engine test::generator = default_random_engine();


void test::testrng()
{
	for (int i = 0; i < 100; i++)
	{
		double r = distr(generator);
		cout << r << endl;
		int s_i_new = (r < 0.5 ? 1 : -1);

		if (r < 0.5) {
			cout << "smaller than 0.5" << endl;
		}
		else {
			cout << "greater than 0.5" << endl;
		}
		
		cout << "s = " << s_i_new << endl;;
	}
	


}
