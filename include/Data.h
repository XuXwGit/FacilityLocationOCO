#ifndef _DATAGENERATOR_H_
#define _DATAGENERATOR_H_

#include "IncludeSTD.h"

class OCOData
{
	typedef vector<vector<double>>										double2;
	typedef vector<vector<vector<double>>>					double3;
	typedef vector<vector<vector<vector<double>>>>	double4;

public:
	OCOData(int num_suppliers, int num_retailers, int num_cross_docks, int num_products, int num_samples)
		: num_suppliers(num_suppliers), num_retailers(num_retailers), num_cross_docks(num_cross_docks), num_products(num_products) , num_samples(num_samples)
	{
		setData();
	}

	void setData();
	void set_demand();
	void set_varation_demand();
	void set_position();
	void set_distance();
	void set_unit_cost();
	void set_trans_time();
	void set_service_level();
	void set_volume();
	void set_truck_capacity();
	void set_C1();
	void set_C2();

	inline int getNumSuppliers() {
		return num_suppliers;
	}
	inline int getNumCrossDocks() {
		return num_cross_docks;
	}
	inline int getNumRetailers() {
		return num_retailers;
	}
	inline int getNumProducts() {
		return num_products;
	}
	inline int getNumSamples() {
		return num_samples;
	}
	inline vector<vector<int>>& getDemand() {
		return demand;
	}
	inline vector<vector<vector<int>>>& getVarationDemand() {
		return varation_demand;
	}
	inline vector<double>& getReservationCost() {
		return unit_reservation_cost;
	}
	inline double4& getC1() {
		return C1;
	}
	inline double4& getC2() {
		return C2;
	}
	inline double2& getBeta() {
		return beta;
	}
	//~OCOData();

private:
	// problem scale:
	// |I|, |J|, |K|, |P|
	int num_suppliers;
	int num_retailers;
	int num_cross_docks;
	int num_products;
	// N =  |I| * |J| * |K| * |P|
	// int N;

	// S
	int num_samples = 100;

	// base demand
	// demand : D[j][p] for j \in J, p \in P
	// |J| x |P|
	vector<vector<int>> demand;
	// demand ~ U[lp, up]
	int demand_lower_bound = 5;
	int demand_upper_bound = 50;

	// Varation Demand
	vector<vector<vector<int>>> varation_demand;
	double uncertain_degree = 0.10;

	// position (longitude, latitude)
	// |I| x 2
	vector<pair<double, double>> position_of_suppliers;
	// |J| X 2
	vector<pair<double, double>> position_of_retailers;
	// |P| X 2
	vector<pair<double, double>> position_of_crossdock;
	// longitude/latitude ~ U[lp, up]
	int longitude_lower_bound = 0;
	int longitude_upper_bound = 1000;
	int latitude_lower_bound = 0;
	int latitude_upper_bound = 1000;

	// distance
	// |I| X |K|
	vector<vector<double>> supplier_to_crossdock;
	// |K| X |J|
	vector<vector<double>> crossdock_to_retailer;


	// unit cost 
	// f_k : unit reservation cost for cross-dock k
	// |K|
	vector<double> unit_reservation_cost;
	// f_k ~ U[lp, up]
	int unit_reservation_cost_lower_bound = 20;
	int unit_reservation_cost_upper_bound = 40;

	// c_ik_I: unit transportation cost from supplier i to cross-dock k
	// |I| X |K|
	vector<vector<double>> unit_transportation_cost_I;
	// c_kj_O: unit transportation cost from cross-dock k to retailer j 
	// |K| X |J|
	vector<vector<double>> unit_transportation_cost_O;
	// h_p : unit inventory holding cost of product p
	// |P|

	// C_ijkp_1
	// C_ijkp_2
	double4 C1;
	double4 C2;

	vector<double> unit_inventory_cost;
	// h_p ~ U[lp, up]
	int unit_inventory_cost_lower_bound = 1;
	int unit_inventory_cost_upper_bound = 10;

	// transportation time
	// distance / 4000
	// |I| X |K|
	vector<vector<double>> trans_time_supplier_to_crossdock;
	// |K| X |J|
	vector<vector<double>> trans_time_crossdock_to_retailer;
	// speed of truck : 4000
	int truck_speed = 4000;

	// service level
	// beta[j][p]
	// |J| X |P|
	vector<vector<double>> beta;
	// beta ~ U[1, 10]
	double beta_lower_bound = 0.8;
	double beta_upper_bound = 1.0;

	// volume of product p : v_p
	// |P|
	vector<int> product_volume;
	// volume ~ U[1, 10]
	int volume_lower_bound = 1;
	int volume_upper_bound = 10;

	  
	// truck capacity
	// V
	int truck_capacity;
	// capacity ~ U[100, 1000]
	int capacity_lower_bound = 100;
	int capacity_upper_bound = 1000;

	// random engine
	default_random_engine e;
};

#endif // !_DATAGENERATOR_H_
