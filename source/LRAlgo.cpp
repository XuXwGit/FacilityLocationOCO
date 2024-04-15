#include "LRAlgo.h"

void LRAlgo::frame(int type)
{
	initialize();

	double Gap = epsilon;
	int n = 0;
	while (Gap >= epsilon) {
		n++;
		if (type == 1) {
			calculate_XY_directly();
		}
		else if (type == 2)
		{
			calculate_XY_byOCO();
		}
		update_lambda(n);
		
		double a = *(max_element(delta.begin(), delta.end()));
		double b = *(min_element(delta.begin(), delta.end()));
		Gap = abs(a) > abs(b) ? abs(a) : abs(b);
		cout << " Gap " << n << " : " << Gap <<"\t" <<  Y_value[0] <<"\t" << current_obj << "\t" << optimal_obj << "\t" << endl;
	}
	cout << endl;

	printsolution();
}

void LRAlgo::initialize()
{
	// initialize lambda
	lambda = vector<vector<double>>(Data.getNumRetailers(), vector<double>(Data.getNumProducts(), double(0.0)));

	// initialize delta = lambda^n+1 - lambda^n
	delta = vector<double>(Data.getNumRetailers() * Data.getNumProducts());
	for (size_t j = 0; j < Data.getNumRetailers(); j++)
	{
		for (size_t p = 0; p < Data.getNumProducts(); p++)
		{
			delta[j * Data.getNumProducts() + p] = lambda[j][p];
		}
	}

	// calculate average demand
	// beta[j][p]
	vector<vector<double>> beta = Data.getBeta();
	// d[j][p]
	vector<vector<int>> base_demand = Data.getDemand();
	// vd[j][p][s]
	vector<vector<vector<int>>> varation_demand = Data.getVarationDemand();
	average_service_demand = vector<vector<double>>(Data.getNumRetailers(), vector<double>(Data.getNumProducts(), 0));
	for (size_t j = 0; j < Data.getNumRetailers(); j++)
	{
		for (size_t p = 0; p < Data.getNumProducts(); p++)
		{
			for (size_t s = 0; s < Data.getNumSamples(); s++)
			{
				average_service_demand[j][p] += beta[j][p] * (base_demand[j][p] + varation_demand[s][j][p]) / Data.getNumSamples();
			}
		}
	}

	current_obj = 0;
	optimal_obj = 0;
}

void LRAlgo::calculate_XY_directly()
{
	LRM.reset_lambda_in_obj(lambda);
	LRM.solve_model();
	Y_value = LRM.get_Yvalue();
	X_value = LRM.get_Xvalue();
	current_obj = LRM.getObjVal();
	if (optimal_obj < current_obj) {
		optimal_obj = current_obj;
	}
}

void LRAlgo::calculate_XY_byOCO()
{
	OCO.frame(lambda);
	Y_value = OCO.getYvalue();
	X_value = OCO.getXvalue();
}

void LRAlgo::update_lambda(int n)
{
	// update delta
	// step = 1 / sqrt(n)
	double step = 1 / sqrt(n);

	delta.clear();
	delta = vector<double>(Data.getNumRetailers() * Data.getNumProducts(), 0);
	for (size_t j = 0; j < Data.getNumRetailers(); j++)
	{
		for (size_t p = 0; p < Data.getNumProducts(); p++)
		{
			delta[j * Data.getNumProducts() + p] = average_service_demand[j][p];
			for (size_t i = 0; i < Data.getNumSuppliers(); i++)
			{
				for (size_t k = 0; k < Data.getNumCrossDocks(); k++)
				{
					delta[j * Data.getNumProducts() + p] -= X_value[i][j][k][p];
				}
			}
			delta[j * Data.getNumProducts() + p] = delta[j * Data.getNumProducts() + p] * step;
		}
	}

	// update lambda
	for (size_t j = 0; j < Data.getNumSuppliers(); j++)
	{
		for (size_t p = 0; p < Data.getNumProducts(); p++)
		{
			lambda[j][p] += delta[j * Data.getNumProducts() + p];
		}
	}

}

void LRAlgo::printsolution()
{
	cout << "Y value : ";
	for (size_t k = 0; k < Data.getNumCrossDocks(); k++)
	{
		cout << Y_value[k] << "\t";
	}
}
