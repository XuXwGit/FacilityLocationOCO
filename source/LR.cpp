#include <LR.h>

LRModel::LRModel(OCOData& data)
	: Data(data)
{
	this->env = IloEnv();
	this->model = IloModel(env);
	this->cplex = IloCplex(model);

	lambda = vector<vector<double>>(Data.getNumRetailers(), vector<double>(Data.getNumProducts(), 0.00));

	createDecisions();
	setObjective();
	addConstraints();
}

void LRModel::createDecisions()
{
	// create Stage I decisions : reservation capacity
	// Y[k] for k \in K
	Y = IloNumVarArray(env, Data.getNumCrossDocks(), 0, IloInfinity, IloNumVar::Float);
	for (size_t k = 0; k < Data.getNumCrossDocks(); k++)
	{
		Y[k].setName(("Y" + to_string(k)).c_str());
	}

	// create Stage II decisions : flow distribution
	// X[i][j][k][p][s]
	X = IloNumVar5(env, Data.getNumSuppliers());
	for (size_t i = 0; i < Data.getNumSuppliers(); i++)
	{
		X[i] = IloNumVar4(env, Data.getNumRetailers());
		for (size_t j = 0; j < Data.getNumRetailers(); j++)
		{
			X[i][j] = IloNumVar3(env, Data.getNumCrossDocks());
			for (size_t k = 0; k < Data.getNumCrossDocks(); k++)
			{
				X[i][j][k] = IloNumVar2(env, Data.getNumProducts());
				for (size_t p = 0; p < Data.getNumProducts(); p++)
				{
					X[i][j][k][p] = IloNumVarArray(env, Data.getNumSamples(), 0, IloInfinity, IloNumVar::Float);
					for (size_t s = 0; s < Data.getNumSamples(); s++)
					{
						X[i][j][k][p][s].setName(("X" + to_string(i) + "_" + to_string(j)
							+ "_" + to_string(p) + "_" + to_string(k)
							+ "_" + to_string(s)).c_str());
					}
				}
			}
		}
	}
}

void LRModel::setObjective()
{
	vector<double> F_k = Data.getReservationCost();
	double4 C1 = Data.getC1();
	double4 C2 = Data.getC2();

	IloExpr obj(env);
	// add item 1
	for (size_t k = 0; k < Data.getNumCrossDocks(); k++)
	{
		obj += F_k[k] * Y[k];
	}

	// add item 2
	for (size_t s = 0; s < Data.getNumSamples(); s++)
	{
		for (size_t i = 0; i < Data.getNumSuppliers(); i++)
		{
			for (size_t j = 0; j < Data.getNumRetailers(); j++)
			{
				for (size_t k = 0; k < Data.getNumCrossDocks(); k++)
				{
					for (size_t p = 0; p < Data.getNumProducts(); p++)
					{
						obj += (C1[i][j][k][p] + C2[i][j][k][p]) / Data.getNumSamples() * X[i][j][k][p][s];
					}
				}
			}
		}
	}

	// beta[j][p]
	vector<vector<double>> beta = Data.getBeta();
	// d[j][p]
	vector<vector<int>> base_demand = Data.getDemand();
	// vd[j][p][s]
	vector<vector<vector<int>>> varation_demand = Data.getVarationDemand();
	for (size_t j = 0; j < Data.getNumRetailers(); j++)
	{
		for (size_t p = 0; p < Data.getNumProducts(); p++)
		{
			for (size_t s = 0; s < Data.getNumSamples(); s++)
			{
				obj += lambda[j][p] * beta[j][p] * (base_demand[j][p] + varation_demand[s][j][p]) / Data.getNumSamples();

				for (size_t i = 0; i < Data.getNumSuppliers(); i++)
				{
					for (size_t k = 0; k < Data.getNumCrossDocks(); k++)
					{
						obj -= lambda[j][p] * X[i][j][k][p][s] / Data.getNumSamples();
					}
				}
			}
		}
	}

	objective = IloMinimize(env, obj);

	model.add(objective);
}

void LRModel::addConstraint2()
{
	// d[j][p]
	vector<vector<int>> base_demand = Data.getDemand();
	// vd[s][j][p]
	vector<vector<vector<int>>> varation_demand = Data.getVarationDemand();

	for (size_t j = 0; j < Data.getNumRetailers(); j++)
	{
		for (size_t p = 0; p < Data.getNumProducts(); p++)
		{
			for (size_t s = 0; s < Data.getNumSamples(); s++)
			{
				IloExpr left(env);

				for (size_t i = 0; i < Data.getNumSuppliers(); i++)
				{
					for (size_t k = 0; k < Data.getNumCrossDocks(); k++)
					{
						left += X[i][j][k][p][s];
					}
				}

				string constr = "C2_" + to_string(j) + "_" + to_string(p) + "_" + to_string(s);
				model.add(left <= base_demand[j][p] + varation_demand[s][j][p]).setName(constr.c_str());
			}
		}
	}
}

void LRModel::addConstraint3()
{
	for (size_t k = 0; k < Data.getNumCrossDocks(); k++)
	{
		for (size_t s = 0; s < Data.getNumSamples(); s++)
		{
			IloExpr left(env);

			for (size_t i = 0; i < Data.getNumSuppliers(); i++)
			{
				for (size_t j = 0; j < Data.getNumRetailers(); j++)
				{
					for (size_t p = 0; p < Data.getNumProducts(); p++)
					{
						left += X[i][j][k][p][s];
					}
				}
			}

			string constr = "C3_" + to_string(k) + "_" + to_string(s);
			model.add(left <= Y[k]).setName(constr.c_str());
		}
	}
}

void LRModel::reset_lambda_in_obj(vector<vector<double>>& lambda_value)
{
	 lambda = lambda_value;

	 vector<double> F_k = Data.getReservationCost();
	 // beta[j][p]
	 vector<vector<double>> beta = Data.getBeta();
	 // d[j][p]
	 vector<vector<int>> base_demand = Data.getDemand();
	 // vd[j][p][s]
	 vector<vector<vector<int>>> varation_demand = Data.getVarationDemand();
	 // C1
	 double4 C1 = Data.getC1();
	 // C2
	 double4 C2 = Data.getC2();

	 // calculate coeff of X[i][j][k][p][s]
	 for (size_t s = 0; s < Data.getNumSamples(); s++)
	 {
		 for (size_t i = 0; i < Data.getNumSuppliers(); i++)
		 {
			 for (size_t j = 0; j < Data.getNumRetailers(); j++)
			 {
				 for (size_t k = 0; k < Data.getNumCrossDocks(); k++)
				 {
					 for (size_t p = 0; p < Data.getNumProducts(); p++)
					 {
						 IloNum new_coeff = (C1[i][j][k][p] + C2[i][j][k][p] - lambda[j][p]) / Data.getNumSamples();
						 objective.setLinearCoef(X[i][j][k][p][s], new_coeff);
					 }
				 }
			 }
		 }
	 }

	 // calculate coeff for constant
	 IloNum constant_coeff = 0;
	 for (size_t j = 0; j < Data.getNumRetailers(); j++)
	 {
		 for (size_t p = 0; p < Data.getNumProducts(); p++)
		 {
			 for (size_t s = 0; s < Data.getNumSamples(); s++)
			 {
				 constant_coeff += lambda[j][p] * beta[j][p] * (base_demand[j][p] + varation_demand[s][j][p]) / Data.getNumSamples();
			 }
		 }
	 }
	 objective.setConstant(constant_coeff);
}

void LRModel::solve_model()
{
	string model_name = "LR" + to_string(Data.getNumSuppliers())
		+ "_" + to_string(Data.getNumRetailers())
		+ "_" + to_string(Data.getNumCrossDocks())
		+ "_" + to_string(Data.getNumProducts())
		+ "_" + to_string(Data.getNumSamples()) + ".lp";
	cplex.exportModel(model_name.c_str());

	if (cplex.solve()) {
		setXvalue();
		setYvalue();
		cplex.setOut(env.getNullStream());
		//printsolution();
	}
}

void LRModel::setYvalue()
{
	Yvalue = vector<double>(Data.getNumCrossDocks(), 0.0);
	for (size_t k = 0; k < Data.getNumCrossDocks(); k++)
	{
		Yvalue[k] = cplex.getValue(Y[k]);
	}
}

void LRModel::setXvalue()
{
	Xvalue.clear();
	Xvalue = double4(Data.getNumSuppliers(), double3(Data.getNumRetailers(),
		double2(Data.getNumCrossDocks(), vector<double>(Data.getNumProducts(),
			double(0)))));
	for (size_t i = 0; i < Data.getNumSuppliers(); i++)
	{
		for (size_t j = 0; j < Data.getNumRetailers(); j++)
		{
			for (size_t k = 0; k < Data.getNumCrossDocks(); k++)
			{
				for (size_t p = 0; p < Data.getNumProducts(); p++)
				{
					for (size_t s = 0; s < Data.getNumSamples(); s++)
					{
						Xvalue[i][j][k][p] += cplex.getValue(X[i][j][k][p][s]) / Data.getNumSamples();
					}
				}
			}
		}
	}
}

void LRModel::printsolution()
{
	cout << "Solution status: " << cplex.getStatus() << endl;
	cout << "The Optimal Objective: " << cplex.getObjValue() << endl;
	cout << "Solving Time: " << cplex.getCplexTime() << endl;
}