#include "OGD.h"

OGD::OGD(OCOData& Data)
	: Data(Data), DSP(DualSubLR(Data))
{
	initialize();
}

OGD::OGD(OCOData& Data, vector<vector<double>>& lambda_value)
	: Data(Data), lambda(lambda_value), DSP(DualSubLR(Data, lambda_value))
{
	initialize();
}

void OGD::initialize()
{
	Yvalue = vector<double>(Data.getNumCrossDocks(), 0.0);
	Xvalue = double4(Data.getNumSuppliers(), double3(Data.getNumRetailers(), double2(Data.getNumCrossDocks(), vector<double>(Data.getNumProducts(), 0.0))));
}

void OGD::frame(vector<vector<double>>& lambda)
{
	update_multiplier(lambda);
	DSP.update_lambda_in_model(lambda);

	std::fill(Yvalue.begin(), Yvalue.end(), 0);
	std::fill(Xvalue.begin(), Xvalue.end(), double3(Data.getNumRetailers(), double2(Data.getNumCrossDocks(), vector<double>(Data.getNumProducts(), 0))));
	for (int s = 0; s < Data.getNumSamples(); s++)
	{
		DSP.update_obj(s);
		DSP.solve_model();
		//Yset.push_back(DSP.getYvalue());
		//Xset.push_back(DSP.getXvalue());

		for (size_t k = 0; k < Data.getNumCrossDocks(); k++)
		{
			Yvalue[k] += DSP.getYvalue()[k] / Data.getNumSamples();
		}

		for (size_t i = 0; i < Data.getNumSuppliers(); i++)
		{
			for (size_t j = 0; j < Data.getNumRetailers(); j++)
			{
				for (size_t k = 0; k < Data.getNumCrossDocks(); k++)
				{
					for (size_t p = 0; p < Data.getNumProducts(); p++)
					{
						Xvalue[i][j][k][p] += DSP.getXvalue()[i][j][k][p] / Data.getNumSamples();
					}
				}
			}
		}
	}


	DSP.setY0(Yvalue);
}

void OGD::update_multiplier(vector<vector<double>>& lambda_value)
{
	lambda = lambda_value;
}

