#include "DualSubLR.h"

DualSubLR::DualSubLR(OCOData& Data)
	: Data(Data), lambda(vector<vector<double>>(Data.getNumRetailers(), vector<double>(Data.getNumProducts(), 0.00)))
{
	initialize();
}

DualSubLR::DualSubLR(OCOData& Data, vector<vector<double>>& lambda)
	: Data(Data), lambda(lambda)
{
	initialize();
}

void DualSubLR::initialize()
{
	this->env = IloEnv();
	this->model = IloModel(env);
	this->cplex = IloCplex(model);

	Yvalue = vector<double>(Data.getNumCrossDocks(), 0.0);
	Xvalue = double4(Data.getNumSuppliers(), double3(Data.getNumRetailers(),
						double2(Data.getNumCrossDocks(),
							vector<double>(Data.getNumProducts(), double(0.0)))));

	createDecisions();
	setObjective();
	addConstraints();
}

void DualSubLR::createDecisions()
{
	// create Stage I decisions : reservation capacity
	// gamma[k] for k \in K
	gamma = IloNumVarArray(env, Data.getNumCrossDocks(), -IloInfinity, 0, IloNumVar::Float);
	for (size_t k = 0; k < Data.getNumCrossDocks(); k++)
	{
		gamma[k].setName(("gamma" + to_string(k)).c_str());
	}

	// create Stage II decisions : flow distribution
	// alpha[j][p]
	alpha = IloArray<IloNumVarArray>(env, Data.getNumRetailers());
	for (size_t j = 0; j < Data.getNumRetailers(); j++)
	{
		alpha[j] = IloNumVarArray(env, Data.getNumProducts());
		for (size_t p = 0; p < Data.getNumProducts(); p++)
		{
			alpha[j][p] = IloNumVar(env, -IloInfinity, 0, IloNumVar::Float);
			alpha[j][p].setName(("alpha" + to_string(j) + "_" + to_string(p)).c_str());
		}
	}
}

void DualSubLR::setObjective()
{
	// d[j][p]
	vector<vector<int>> base_demand = Data.getDemand();
	// vd[j][p][s]
	vector<vector<vector<int>>> varation_demand = Data.getVarationDemand();

	IloExpr obj(env);

	// add item 1
	for (size_t j = 0; j < Data.getNumRetailers(); j++)
	{
		for (size_t p = 0; p < Data.getNumProducts(); p++)
		{
			obj += (base_demand[j][p] + varation_demand[scene_index][j][p]) * alpha[j][p];
		}
	}
	
	// add item 2
	for (size_t k = 0; k < Data.getNumCrossDocks(); k++)
	{
		obj += Yvalue[k] * gamma[k];
	}

	objective = IloMaximize(env, obj);
	model.add(objective);
}

void DualSubLR::addConstraints()
{
	double4 C1 = Data.getC1();
	double4 C2 = Data.getC2();

	constraints = range4(Data.getNumSuppliers(), range3(Data.getNumRetailers(), range2(Data.getNumCrossDocks(), vector<IloRange>(Data.getNumProducts(), IloRange()))));
	for (size_t i = 0; i < Data.getNumSuppliers(); i++)
	{
		for (size_t j = 0; j < Data.getNumRetailers(); j++)
		{
			for (size_t k = 0; k < Data.getNumCrossDocks(); k++)
			{
				for (size_t p = 0; p < Data.getNumProducts(); p++)
				{
					constraints[i][j][k][p] = IloRange(env, alpha[j][p] + gamma[k] , C1[i][j][k][p] + C2[i][j][k][p] - lambda[j][p]);
					model.add(constraints[i][j][k][p]);
				}
			}
		}
	}
}

void DualSubLR::update_lambda_in_model(vector<vector<double>>& lambda_value)
{
	lambda = lambda_value;

	double4 C1 = Data.getC1();
	double4 C2 = Data.getC2();

	for (size_t i = 0; i < Data.getNumSuppliers(); i++)
	{
		for (size_t j = 0; j < Data.getNumRetailers(); j++)
		{
			for (size_t k = 0; k < Data.getNumCrossDocks(); k++)
			{
				for (size_t p = 0; p < Data.getNumProducts(); p++)
				{
					constraints[i][j][k][p].setUB(C1[i][j][k][p] + C2[i][j][k][p] - lambda[j][p]);
				}
			}
		}
	}
}

void DualSubLR::update_obj(int scene_index)
{
	set_scene_index(scene_index);
	
	// d[j][p]
	vector<vector<int>> base_demand = Data.getDemand();
	// vd[j][p][s]
	vector<vector<vector<int>>> varation_demand = Data.getVarationDemand();

	// add item 1
	for (size_t j = 0; j < Data.getNumRetailers(); j++)
	{
		for (size_t p = 0; p < Data.getNumProducts(); p++)
		{
			IloNum new_coeff = base_demand[j][p] + varation_demand[scene_index][j][p];
			objective.setLinearCoef(alpha[j][p], new_coeff);
		}
	}

	// add item 2
	for (size_t k = 0; k < Data.getNumCrossDocks(); k++)
	{
		objective.setLinearCoef(gamma[k], Yvalue[k]);
	}
}

void DualSubLR::solve_model()
{
	cplex.setOut(env.getNullStream());
	if (cplex.solve()) {
		updateYvalue();
		setXvalue();
		//printsolution();
	}
	else
	{
		cplex.exportModel("error.lp");
		cout << cplex.getStatus() << endl;
		cout << "There exits an error!" << endl;
	}
}

// update Yvalue = max{Y_k^n(s-1) - ¦Ç(s) (f_k + r_k(s)), 0}
void DualSubLR::updateYvalue()
{
	vector<double> F_k = Data.getReservationCost();

	for (size_t k = 0; k < Data.getNumCrossDocks(); k++)
	{
		double tempYvalue = Yvalue[k] -(F_k[k] + cplex.getValue(gamma[k])) / sqrt(scene_index + 1);
		Yvalue[k] =  tempYvalue > 0 ? tempYvalue : 0;
	}
}

void DualSubLR::setXvalue()
{
	for (size_t i = 0; i < Data.getNumSuppliers(); i++)
	{
		for (size_t j = 0; j < Data.getNumRetailers(); j++)
		{
			for (size_t k = 0; k < Data.getNumCrossDocks(); k++)
			{
				for (size_t p = 0; p < Data.getNumProducts(); p++)
				{
					Xvalue[i][j][k][p] = cplex.getDual(constraints[i][j][k][p]);
				}
			}
		}
	}
}

void DualSubLR::printsolution()
{
	cout << "Solution status: " << cplex.getStatus() << endl;
	cout << "The Optimal Objective: " << cplex.getObjValue() << endl;
	cout << "Solving Time: " << cplex.getCplexTime() << endl;
}
