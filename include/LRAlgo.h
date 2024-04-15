#pragma once
#ifndef _LR_ALGO_H_
#define _LR_ALGO_H_

#include "IncludeSTD.h"
#include "Data.h"
#include "LR.h"
#include <OGD.h>

class LRAlgo
{
	typedef vector<vector<vector<vector<double>>>> double4;

public:
	inline LRAlgo(OCOData& Data, int type) 
		: Data(Data), 
		LRM(LRModel(Data)), 
		OCO(OGD(Data))
	{
		frame(type);
	}
	
	void frame(int type);

	void initialize();
	void calculate_XY_directly();
	void calculate_XY_byOCO();
	void update_lambda(int n);

	void printsolution();

private:
	// Lagrangian dual multiplier \lambda_{jp}
	vector<vector<double>> lambda;
	vector<double> Y_value;
	double4 X_value;

	OCOData Data;
	LRModel LRM;
	OGD OCO;

	vector<vector<double>> average_service_demand;

	// gap = lambda^(n+1) - lambda^n
	vector<double> delta;
	double epsilon = 0.01;
	double current_obj;
	double optimal_obj;
public:

	LRAlgo() = default;
};

#endif // !_LR_ALGO_H_
