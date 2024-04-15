#pragma once
#ifndef _LR_H_
#define _LR_H_

#include "IncludeSTD.h"

#include <Data.h>

class LRModel
{
	typedef IloArray<IloArray<IloArray<IloArray<IloNumVarArray>>>>	IloNumVar5;
	typedef IloArray<IloArray<IloArray<IloNumVarArray>>>						IloNumVar4;
	typedef IloArray<IloArray<IloNumVarArray>>												IloNumVar3;
	typedef IloArray<IloNumVarArray>																	IloNumVar2;
	typedef vector<vector<vector<vector<double>>>>									double4;
	typedef vector<vector<vector<double>>>													double3;
	typedef vector<vector<double>>																		double2;
public:
	LRModel(OCOData& data);

	void createDecisions();
	void setObjective();
	inline void addConstraints() {
		addConstraint2();
		addConstraint3();
	}
	// demand constraint
	void addConstraint2();
	// capacity constraint
	void addConstraint3();

	void reset_lambda_in_obj(vector<vector<double>>& lambda_value);
	void solve_model();

	void setYvalue();
	void setXvalue();

	inline vector<double>& get_Yvalue() {
		return Yvalue;
	}
	inline double4& get_Xvalue() {
		return Xvalue;
	}

	inline double getObjVal() {
		return cplex.getObjValue();
	}

	void printsolution();

private:
	// data
	OCOData Data;

	// cplex model
	IloEnv env;
	IloModel model;
	IloCplex cplex;
	IloObjective objective;

	IloExpr fixed_item_in_obj;

	// decision variables
		// Y_k for k \in K
	IloNumVarArray Y;
	// X_ijkp^s for i \in I, j \in J, k \in K, p \in P for s \in S
	IloArray<IloArray<IloArray<IloArray<IloNumVarArray>>>> X;

	vector<double> Yvalue;
	double4	 Xvalue;

	// multiplier
	vector<vector<double>> lambda;
public:

	LRModel() = default;
};
#endif // !_LR_H_
