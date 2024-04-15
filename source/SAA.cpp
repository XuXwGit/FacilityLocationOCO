#include "SAA.h"

SAAModel::SAAModel(OCOData& data)
	: Data(data), num_samples(Data.getNumSamples())
{
	this->env = IloEnv();
	this->model = IloModel(env);
	this->cplex = IloCplex(model);

	createDecisions();
	setObjective();
	addConstraints();
	solveModel();
}

void SAAModel::createDecisions()
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
					X[i][j][k][p] = IloNumVarArray(env, num_samples, 0, IloInfinity, IloNumVar::Float);
					for (size_t s = 0; s < num_samples; s++)
					{
						X[i][j][k][p][s].setName(("X" + to_string(i) + "_" + to_string(j) 
																		+ "_" + to_string(p) + "_" + to_string(k) 
																		+ "_" + to_string(s)										).c_str());
					}
				}
			}
		}
	}
}

void SAAModel::setObjective()
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
	for (size_t s = 0; s < num_samples; s++)
	{
		for (size_t i = 0; i < Data.getNumSuppliers(); i++)
		{
			for (size_t j = 0; j < Data.getNumRetailers(); j++)
			{
				for (size_t k = 0; k < Data.getNumCrossDocks(); k++)
				{
					for (size_t p = 0; p < Data.getNumProducts(); p++)
					{
						obj += (C1[i][j][k][p] + C2[i][j][k][p]) / num_samples * X[i][j][k][p][s];
					}
				}
			}
		}
	}

	model.add(IloMinimize(env, obj));
}

void SAAModel::addConstraint1()
{
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
			IloExpr left(env);
			for (size_t s = 0; s < num_samples; s++)
			{
				for (size_t i = 0; i < Data.getNumSuppliers(); i++)
				{
					for (size_t k = 0; k < Data.getNumCrossDocks(); k++)
					{
						left += X[i][j][k][p][s];
					}
				}

				left += - beta[j][p] * (base_demand[j][p] + varation_demand[s][j][p]);
			}

			string constr = "C1_" + to_string(j) + "_" + to_string(p);
			model.add(left >= 0).setName(constr.c_str());
		}
	}
}

void SAAModel::addConstraint2()
{
	// d[j][p]
	vector<vector<int>> base_demand = Data.getDemand();
	// vd[s][j][p]
	vector<vector<vector<int>>> varation_demand = Data.getVarationDemand();

	for (size_t j = 0; j < Data.getNumRetailers(); j++)
	{
		for (size_t p = 0; p < Data.getNumProducts(); p++)
		{
			for (size_t s = 0; s < num_samples; s++)
			{
				IloExpr left(env);

				for (size_t i = 0; i < Data.getNumSuppliers(); i++)
				{
					for (size_t k = 0; k < Data.getNumCrossDocks(); k++)
					{
						left += X[i][j][k][p][s];
					}
				}

				string constr = "C2_" + to_string(j) + "_" + to_string(p)+"_"+to_string(s);
				model.add(left <= base_demand[j][p] + varation_demand[s][j][p]).setName(constr.c_str());
			}
		}
	}
}

void SAAModel::addConstraint3()
{
	for (size_t k = 0; k < Data.getNumCrossDocks(); k++)
	{
		for (size_t s = 0; s < num_samples; s++)
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

void SAAModel::solveModel()
{
	string model_name = "SAA" + to_string(Data.getNumSuppliers()) 
												+ "_" + to_string(Data.getNumRetailers())
												+ "_" + to_string(Data.getNumCrossDocks())
												+ "_" + to_string(Data.getNumProducts())
												+ "_" + to_string(num_samples) + ".lp";
	cplex.exportModel(model_name.c_str());
	cplex.setOut(env.getNullStream());
	if (cplex.solve()) {
		Y_value = vector<double>(Data.getNumCrossDocks(), 0);
		for (size_t k = 0; k < Data.getNumCrossDocks(); k++)
		{
			Y_value[k] = cplex.getValue(Y[k]);
		}

		printsolution();
	}
}

void SAAModel::printsolution()
{
	cout << "Solution status: " << cplex.getStatus() << endl;
	cout << "The Optimal Objective: " << cplex.getObjValue() << endl;
	cout << "Solving Time: " << cplex.getCplexTime() << endl;

	for (size_t k = 0; k < Data.getNumCrossDocks(); k++)
	{
		cout <<"Y value : " << Y_value[k] << "\t";
	}
	cout << endl;
}

