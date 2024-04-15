#pragma once
#ifndef _SAA_H_
#define _SAA_H_

#include "Data.h"

/// <summary>
/// This class can bulid a stochastic model 
/// with sample approximate average (SAA) method.
/// </summary>
class SAAModel
{
	typedef IloArray<IloArray<IloArray<IloArray<IloNumVarArray>>>>	IloNumVar5;
	typedef IloArray<IloArray<IloArray<IloNumVarArray>>>						IloNumVar4;
	typedef IloArray<IloArray<IloNumVarArray>>												IloNumVar3;
	typedef IloArray<IloNumVarArray>																	IloNumVar2;

	typedef vector<vector<vector<vector<double>>>>									double4;

public:
	SAAModel(OCOData& data);

	void createDecisions();
	void setObjective();
	inline void addConstraints() {
		addConstraint1();
		addConstraint2();
		addConstraint3();
	}

	// service level constraint
	void addConstraint1();
	// demand constraint
	void addConstraint2();
	// capacity constraint
	void addConstraint3();

	void solveModel();

	void printsolution();

private:
	OCOData Data;
	// S
	int num_samples =100;

	// cplex model
	IloEnv env;
	IloModel model;
	IloCplex cplex;
	// Y_k for k \in K
	IloNumVarArray Y;
	// X_ijkp^s for i \in I, j \in J, k \in K, p \in P for s \in S
	IloArray<IloArray<IloArray<IloArray<IloNumVarArray>>>> X;

	// solution
	vector<double> Y_value;
};

#endif // !_SAA_H_
