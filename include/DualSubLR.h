#pragma once
#ifndef _DUAL_SUB_LR_H_
#define _DUAL_SUB_LR_H_

#include <Data.h>

class DualSubLR
{
	typedef vector<vector<vector<vector<double>>>>			double4;
	typedef vector<vector<vector<double>>>			double3;
	typedef vector<vector<double>>			double2;

	typedef vector<vector<vector<vector<IloRange>>>>		range4;
	typedef vector<vector<vector<IloRange>>>						range3;
	typedef vector<vector<IloRange>>											range2;
public:
	DualSubLR(OCOData& Data);
	DualSubLR(OCOData& Data, vector<vector<double>>& lambda);
	
	void initialize();
	void createDecisions();
	void setObjective();
	void addConstraints();

	void update_lambda_in_model(vector<vector<double>>& lambda);
	void update_obj(int scene_index);
	void solve_model();

	inline void set_scene_index(int s) {
		scene_index = s;
	}
	inline void setY0(vector<double>& Y0) {
		this->Yvalue = Y0;
	}
	void updateYvalue();
	void setXvalue();

	inline vector<double>& getYvalue() {
		return Yvalue;
	}
	inline double4& getXvalue() {
		return Xvalue;
	}

	void printsolution();

private:
	OCOData Data;
	int scene_index = 0;
	vector<vector<double>> lambda;
	
	// cplex model
	IloEnv env;
	IloModel model;
	IloCplex cplex;	
	// objective
	IloObjective objective;
	// decision variables
	IloArray < IloNumVarArray> alpha;
	IloNumVarArray						   gamma;
	// constraint
	range4 constraints;

	// solution
	vector<double> Yvalue;
	double4 Xvalue;
};

#endif // !_DUAL_SUB_LR_H_
