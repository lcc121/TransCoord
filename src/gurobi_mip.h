#ifndef __GUROBI_MIP_H__
#define __GUROBI_MIP_H__

#include "gurobi_c++.h"
#include "config.h"

#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>

typedef vector<vector<bool>> VVB;

class gurobi
{
public:
	gurobi();
	gurobi(VVB &b, vector<double> &e, VVB &hs_support, vector<int> &unique_cnts);
	~gurobi();

public:
	VVB &B;
	vector<double> &ew;
	VVB &hss;
	vector<int> &ucnts;
	const int tn;             // transcript num
	const int en;             // splice graph edge num
	const int hsn;            // hyper-sets support num

	vector<double> ts;        // transcript abundance 
	vector<int> xs;           // transcript kept or not
	vector<double> obj_values;// each round lp obj values
	int score;                // lp result score: -1=no-solution 0:bad-solution 1:tl-sol 2:opt
	string gid;        

private:
	double ALPHA;
	double BETA;
	double coef;

public:
	int lp();


private:
	int lp_min_max_edge_dev(string glog_file, string glp_file, int lp_idx);
	int lp_min_x_sum(string glog_file, string glp_file, double max_dev, int lp_idx);
	int lp_min_edge_dev_abs_sum(string glog_file, string glp_file, int lp_idx);
	int lp_min_edge_dev_quad_sum(string glog_file, string glp_file, int lp_idx);
	int lp_hier_multiobj_dev_x(string glog_file, string glp_file, int lp_idx);
	int lp_blend_multiobj_dev_x(string glog_file, string glp_file, int lp_idx);
	int lp_hier_multiobj_abs_x(string glog_file, string glp_file, float abs_tol, int lp_idx);

	int initial_model(GRBModel &model);
	GRBVar* add_variables(GRBModel &model, int num, char type, string name_prefix);
	double* get_coefs_x(double exponent);
	
	int add_constr_xt(GRBModel &model, GRBVar *vars_t, GRBVar *vars_x);
	int add_constr_edge_dev_ub(GRBModel &model, GRBVar *vars_t, GRBVar &var_s);
	int add_constr_edge_dev_ub(GRBModel &model, GRBVar *vars_t, double max_dev);
	int add_constr_each_edge_dev(GRBModel & model, GRBVar * vars_t, GRBVar * vars_dev);
	int add_constr_edge_dev_abs(GRBModel & model, GRBVar * vars_dev, GRBVar *vars_abs);
	int add_constr_hss(GRBModel &model, GRBVar *vars_t);
	
	int set_obj_min_edge_diff(GRBModel &model, GRBVar &var_s, string glp_file);
	int set_obj_min_x_sum(GRBModel &model, GRBVar *vars_x, string glp_file);
    int set_obj_min_edge_dev_abs_sum(GRBModel &model, GRBVar *vars_abs, string glp_file);
	int set_obj_min_edge_dev_quad_sum(GRBModel &model, GRBVar *vars_dev, string glp_file);
	int set_hier_multiobj_dev_x(GRBModel &model, GRBVar *vars_x, GRBVar &var_s, string glp_file);
    int set_hier_multiobj_abs_x(GRBModel &model, GRBVar *vars_x, GRBVar *vars_abs, float reltol, string glp_file);
	
	int set_rtn_vars_obj(GRBModel & model, GRBVar * vars_t, GRBVar * vars_x, int & rtn);
	int get_overall_stats(GRBModel &model, int lp_idx);
	int set_score(GRBModel & model, int rtn);

	// output result
	int write(ostream & out, GRBModel & model, int rtn);
	int write_lp_status(ostream &out, GRBModel &model);
	int write_t_x_obj(ostream & out, GRBModel &model);

};
#endif // !__GUROBI_MIP_H__
