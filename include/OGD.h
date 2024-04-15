#pragma once
#include <Data.h>
#include <DualSubLR.h>
#ifndef _OGD_H_
#define _OGD_H_

class OGD
{
	typedef vector<vector<vector<vector<vector<double>>>>>					double5;
	typedef vector<vector<vector<vector<double>>>>									double4;
	typedef vector<vector<vector<double>>>													double3;
	typedef vector<vector<double>>																		double2;

public:
	OGD(OCOData& Data);
	OGD(OCOData& Data, vector<vector<double>>& lambda);

	void initialize();
	void frame(vector<vector<double>>& lambda);
	void update_multiplier(vector<vector<double>>& lambda);

	inline vector<double>& getYvalue() {
		return Yvalue;
	}
	inline double4& getXvalue() {
		return Xvalue;
	}

private:
	OCOData Data;

	DualSubLR DSP;

	vector<double> Yvalue;
	vector<vector<double>> Yset;

	vector<vector<double>> lambda;

	double4 Xvalue;
	double5 Xset;
public:

	OGD() = default;
};


#endif // !_OGD_H_
