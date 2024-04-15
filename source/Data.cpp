#include "Data.h"

void OCOData::setData()
{
	set_demand();
	set_varation_demand();
	set_position();
	set_distance();
	set_unit_cost();
	set_trans_time();
	set_service_level();
	set_volume();
	set_truck_capacity();
	set_C1();
	set_C2();
}

void OCOData::set_demand()
{
	vector<vector<int>> demand;
	for (size_t j = 0; j < num_retailers; j++)
	{
		vector<int> subdemand;
		for (size_t p = 0; p < num_products; p++)
		{
			int a = this->demand_lower_bound;
			int b = this->demand_upper_bound;
			subdemand.push_back((rand() % (b - a + 1)) + a);
		}
		demand.push_back(subdemand);
	}

	this->demand = demand;
}

void OCOData::set_varation_demand()
{
	uniform_int_distribution<int> u( demand_lower_bound * uncertain_degree, demand_upper_bound * uncertain_degree);

	// varation_demand [s][j][p]
	vector<vector<vector<int>>> varation_demand_set = vector<vector<vector<int>>>(num_samples, vector<vector<int>>(num_retailers, vector<int>(num_products, 0)));

	for (size_t s = 0; s < num_samples; s++)
	{

		for (size_t j = 0; j < num_retailers; j++)
		{
			for (size_t p = 0; p < num_products; p++)
			{
				varation_demand_set[s][j][p] = u(e);
			}
		}
	}

	this->varation_demand = varation_demand_set;

}

void OCOData::set_position()
{
	// get  a uniform distribution for latitude
	double lat_lb = this->latitude_lower_bound;
	double lat_ub = this->latitude_upper_bound;
	uniform_real_distribution<double> lat_u(lat_lb, lat_ub);
	double long_lb = this->longitude_lower_bound;
	double long_ub = this->longitude_upper_bound;
	uniform_real_distribution<double> long_u(long_lb, long_ub);

	// set position
	double longitude = 0;
	double latitude = 0;

	// set suppilers' position
	vector<pair<double, double>> suppliers_position;
	for (size_t i = 0; i < num_suppliers; i++)
	{
		// position (longitude, latitude)
		longitude = long_u(e);
		latitude = lat_u(e);
		suppliers_position.push_back(make_pair(longitude, latitude));
	}
	this->position_of_suppliers = suppliers_position;


	// set retailers' position
	vector<pair<double, double>> retailers_position;
	for (size_t i = 0; i < num_retailers; i++)
	{
		// position (longitude, latitude)
		longitude = long_u(e);
		latitude = lat_u(e);
		retailers_position.push_back(make_pair(longitude, latitude));
	}
	this->position_of_retailers = retailers_position;


	// set cross-docks' position
	vector<pair<double, double>> crossdocks_position;
	for (size_t i = 0; i < num_cross_docks; i++)
	{
		// position (longitude, latitude)
		longitude = long_u(e);
		latitude = lat_u(e);
		crossdocks_position.push_back(make_pair(longitude, latitude));
	}
	this->position_of_crossdock = crossdocks_position;
}

void OCOData::set_distance()
{
	// set the distance from supplier i to crossdock k for i \in I, k \in K
	vector<vector<double>> supplier_to_crossdock;
	for (size_t i = 0; i < num_suppliers; i++)
	{
		vector<double> distance_i;
		for (size_t k = 0; k < num_cross_docks; k++)
		{
			pair<double, double> position_i = this->position_of_suppliers[i];
			pair<double, double> position_k = this->position_of_crossdock[k];

			// Manhattan distance
			distance_i.push_back(abs(position_i.first - position_k.first) + abs(position_i.second - position_k.second));
		}
		supplier_to_crossdock.push_back(distance_i);
	}
	this->supplier_to_crossdock = supplier_to_crossdock;

	// set the distance from  crossdock k to retailer j for k \in K, j \in J
	vector<vector<double>> crossdock_to_retailer;
	for (size_t k = 0; k < num_cross_docks; k++)
	{
		vector<double> distance_k;
		for (size_t j = 0; j < num_retailers; j++)
		{
			pair<double, double> position_k = this->position_of_crossdock[k];
			pair<double, double> position_j = this->position_of_retailers[j];

			// Manhattan distance
			distance_k.push_back(abs(position_k.first - position_j.first) + abs(position_k.second - position_j.second));
		}
		crossdock_to_retailer.push_back(distance_k);
	}
	this->crossdock_to_retailer = crossdock_to_retailer;
}

void OCOData::set_unit_cost()
{
	// set unit reservasation cost
	// get  a uniform distribution for reservasation cost
	double rc_lb = this->unit_reservation_cost_lower_bound;
	double rc_ub = this->unit_reservation_cost_upper_bound;
	uniform_real_distribution<double> rc_u(rc_lb, rc_ub);
	vector<double> unit_reservation_cost_set;
	for (size_t k = 0; k < num_cross_docks; k++)
	{
		unit_reservation_cost_set.push_back(rc_u(e));
	}
	this->unit_reservation_cost = unit_reservation_cost_set;


	// set transportation cost
	this->unit_transportation_cost_I = this->supplier_to_crossdock;
	this->unit_transportation_cost_O = this->crossdock_to_retailer;


	// set inventory holding cost
	// get  a uniform distribution for inventory holding cost
	double hc_lb = this->unit_inventory_cost_lower_bound;
	double hc_ub = this->unit_inventory_cost_upper_bound;
	uniform_real_distribution<double> hc_u(hc_lb, hc_ub);
	vector<double> unit_inventory_cost_set;
	for (size_t p = 0; p < num_products; p++)
	{
		unit_inventory_cost_set.push_back(hc_u(e));
	}
	this->unit_inventory_cost = unit_inventory_cost_set;
}

void OCOData::set_trans_time()
{
	// set the transportation time from supplier i to crossdock k for i \in I, k \in K
	vector<vector<double>> time_supplier_to_crossdock = this->supplier_to_crossdock;
	for (size_t i = 0; i < num_suppliers; i++)
	{
		for (size_t k = 0; k < num_cross_docks; k++)
		{
			time_supplier_to_crossdock[i][k] = time_supplier_to_crossdock[i][k] / 4000;
		}
	}
	this->trans_time_supplier_to_crossdock = time_supplier_to_crossdock;

	// set the transportation time from  crossdock k to retailer j for k \in K, j \in J
	vector<vector<double>> time_crossdock_to_retailer = this->crossdock_to_retailer;
	for (size_t k = 0; k < num_cross_docks; k++)
	{
		for (size_t j = 0; j < num_retailers; j++)
		{
			time_crossdock_to_retailer[k][j] = time_crossdock_to_retailer[k][j] / 4000;
		}
	}
	this->trans_time_crossdock_to_retailer = time_crossdock_to_retailer;
}

void OCOData::set_service_level()
{
	double beta_lb = this->beta_lower_bound;
	double beta_ub = this->beta_upper_bound;
	uniform_real_distribution<double> u(beta_lb, beta_ub);

	vector<vector<double>> betaset;
	for (size_t j = 0; j < num_retailers; j++)
	{
		vector<double> products;
		for (size_t p = 0; p < num_products; p++)
		{

			products.push_back(u(e));
		}
		betaset.push_back(products);
	}
	this->beta = betaset;
}

void OCOData::set_volume()
{
	int v_lb = this->volume_lower_bound;
	int v_ub = this->volume_upper_bound;
	uniform_real_distribution<double> u(v_lb, v_ub);

	vector<int> volumns;
	for (size_t p = 0; p < num_products; p++)
	{
		volumns.push_back(u(e));
	}
	this->product_volume = volumns;
}

void OCOData::set_truck_capacity()
{
	// [a, b) : (rand() % (b-a))+ a
	// [a, b] : (rand() % (b-a+1))+ a
	// (a, b] : (rand() % (b-a))+ a + 1
	int a = this->capacity_lower_bound;
	int b = this->capacity_upper_bound;
	uniform_int_distribution<int> u(a, b);
	this->truck_capacity = u(e);
}

void OCOData::set_C1()
{
	vector<vector<double>> c_ik_I = this->unit_transportation_cost_I;
	vector<vector<double>> c_kj_O = this->unit_transportation_cost_O;
	vector<int> v_p = this->product_volume;
	int V = this->truck_capacity;

	this->C1 = double4(num_suppliers, double3(num_retailers, double2(num_cross_docks, vector<double>(num_products, double()))));
	for (size_t i = 0; i < num_suppliers; i++)
	{
		for (size_t j = 0; j < num_retailers; j++)
		{
			for (size_t k = 0; k < num_cross_docks; k++)
			{
				for (size_t p = 0; p < num_products; p++)
				{
					C1[i][j][k][p] = (c_ik_I[i][k] + c_kj_O[k][j]) * v_p[p] / V;
				}
			}
		}
	}
}

void OCOData::set_C2()
{
	vector<vector<double>> t_ik_I = this->trans_time_supplier_to_crossdock;
	vector<vector<double>> t_kj_O = this->trans_time_crossdock_to_retailer;
	vector<int> v_p = this->product_volume;
	vector<double> h_p = unit_inventory_cost;
	int V = this->truck_capacity;

	this->C2 = double4(num_suppliers, double3(num_retailers, double2(num_cross_docks, vector<double>(num_products, double()))));
	for (size_t i = 0; i < num_suppliers; i++)
	{
		for (size_t j = 0; j < num_retailers; j++)
		{
			for (size_t k = 0; k < num_cross_docks; k++)
			{
				for (size_t p = 0; p < num_products; p++)
				{
					C2[i][j][k][p] = h_p[p] * v_p[p] * (t_ik_I[i][k] + t_kj_O[k][j]);
				}
			}
		}
	}

}