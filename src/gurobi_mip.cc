#include "gurobi_mip.h"

#include <sstream>
#include <iostream>
#include <algorithm>

using namespace std;

double x_coef_exp = 3;

gurobi::gurobi(VVB &b, vector<double> &e, VVB &hs_support, vector<int> &unique_cnts)
	:B(b), ew(e), hss(hs_support), ucnts(unique_cnts), tn(B.size()), en(B[0].size()), hsn(hss.size())
{
	assert(tn > 0);
	assert(ew.size() == en);
	assert(ucnts.size() == hsn);
	ts.reserve(tn);
	xs.reserve(tn);
	if (tlb == 0.0025)
	{
		ALPHA = 400;
		BETA = 1e7;
		coef = 4;
	}
	else if (tlb == 0.5)
	{
		ALPHA = 2;
		BETA = 1e9;
		coef = 2;
	}
	else if (tlb == 1)
	{
		ALPHA = 1;
		BETA = 1e9;
		coef = 1;
	}
}

gurobi::~gurobi()
{}

int gurobi::lp()
{
	int rtn = -1; //0:got solution, -1:no solution
	int lp1_rtn = -1;
	int lp2_rtn = -1;
	float abs_tol = 0;
	switch (lp_option)
	{
	case 1:
		// lp12: lp_min_max_edge_dev, lp_min_x_sum
		lp1_rtn = lp_min_max_edge_dev(glog_file1, glp_file1, 1);
		if (lp1_rtn == 0) lp1_rtn = lp_min_x_sum(glog_file2, glp_file2, obj_values[0], 2);
		break;
	case 2:
		// lp_min_edge_dev_abs_sum
		lp1_rtn = lp_min_edge_dev_abs_sum(glog_file1, glp_file1, 1);
		break;
	case 3:
		// lp_min_edge_dev_quad_sum
		lp1_rtn = lp_min_edge_dev_quad_sum(glog_file1, glp_file1, 1);
		break;
	case 4:
		// lp_hier_multiobj_dev_x
		lp1_rtn = lp_hier_multiobj_dev_x(glog_file1, glp_file1, 1);
		break;
	case 5:
		// lp_blend_multiobj_dev_x
		lp1_rtn = lp_blend_multiobj_dev_x(glog_file1, glp_file1, 1);
		break;
	case 6:
		// lp_hier_multiobj_abs_x
		if (abs_tol_option == 2) abs_tol = (float)*max_element(ew.begin(), ew.end());
		lp1_rtn = lp_hier_multiobj_abs_x(glog_file1, glp_file1, abs_tol, 1);
		break;
	default:
		break;
	}
	if (lp1_rtn == 0 || lp2_rtn == 0) rtn = 0;
	return rtn;
}

int gurobi::lp_min_max_edge_dev(string glog_file, string glp_file, int lp_idx)
{
	// minimize max-cov-deviation (between the cov of an edge and the cov-sum of transcripts go through it).
	int rtn = -1;// whether got a solution
	try
	{
		GRBEnv *env = new GRBEnv(glog_file);
		GRBModel model = GRBModel(*env);
		model.set(GRB_StringAttr_ModelName, "min_max_edge_dev");
		initial_model(model);

		// add varibles
		GRBVar *vars_t = add_variables(model, tn, GRB_CONTINUOUS, "t");
		GRBVar *vars_x = add_variables(model, tn, GRB_BINARY, "x");
		GRBVar var_s = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "s");

		// add constraints
		add_constr_xt(model, vars_t, vars_x);
		add_constr_edge_dev_ub(model, vars_t, var_s);
		add_constr_hss(model, vars_t);

		// set objective and run model
		set_obj_min_edge_diff(model, var_s, glp_file);

		// get solution status
		set_rtn_vars_obj(model, vars_t, vars_x, rtn);
		get_overall_stats(model, lp_idx);
		set_score(model, rtn);

		// output
		if (verbose >= 2) write(cout, model, rtn);
		if (verbose == -1)
		{
			ofstream fout;
			fout.open(log_file.c_str(), ios_base::app);
			if (fout.fail()) exit(2);
			write(fout, model, rtn);
		}
		
		// Free environment/vars
		delete[] vars_t;
		delete[] vars_x;
		delete env;
	}
	catch (GRBException e)
	{
		cout << "Error code=" << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...)
	{
		cout << "Exception during optimization" << endl;
	}

	return rtn;// 1:have solution, 0:no solution
}

int gurobi::lp_min_x_sum(string glog_file, string glp_file, double max_dev, int lp_idx)
{
	int rtn = -1; // whether got a solution
	try
	{
		GRBEnv *env = new GRBEnv(glog_file);
		GRBModel model = GRBModel(*env);
		model.set(GRB_StringAttr_ModelName, "min_x_sum");
		initial_model(model);

		// add varibles
		GRBVar *vars_t = add_variables(model, tn, GRB_CONTINUOUS, "t");
		GRBVar *vars_x = add_variables(model, tn, GRB_BINARY, "x");

		// add constraints
		add_constr_xt(model, vars_t, vars_x);
		add_constr_edge_dev_ub(model, vars_t, max_dev);
		add_constr_hss(model, vars_t);

		// set objective and run model
		set_obj_min_x_sum(model, vars_x, glp_file);

		// get solution status and save variables & objective value
		set_rtn_vars_obj(model, vars_t, vars_x, rtn);
		get_overall_stats(model, lp_idx);
		set_score(model, rtn);

		// output
		if (verbose >= 2) write(cout, model, rtn);
		if (verbose == -1)
		{
			ofstream fout;
			fout.open(log_file.c_str(), ios_base::app);
			if (fout.fail()) exit(2);
			write(fout, model, rtn);
		}

		// Free environment/vars
		delete[] vars_t;
		delete[] vars_x;
		delete env;
	}
	catch (GRBException e)
	{
		cout << "Error code=" << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...)
	{
		cout << "Exception during optimization" << endl;
	}
	
	return rtn;
}

int gurobi::lp_min_edge_dev_abs_sum(string glog_file, string glp_file, int lp_idx)
{
	// minimize edge-cov-deviation-abs-sum (between the cov of an edge and the cov-sum of transcripts go through it).
	int rtn = -1; // whether got a solution
	try
	{
		GRBEnv *env = new GRBEnv(glog_file);
		GRBModel model = GRBModel(*env);
		model.set(GRB_StringAttr_ModelName, "min_edge_dev_abs_sum");
		initial_model(model);

		// add varibles
		GRBVar *vars_t = add_variables(model, tn, GRB_CONTINUOUS, "t");
		GRBVar *vars_x = add_variables(model, tn, GRB_BINARY, "x");
		GRBVar *vars_dev = add_variables(model, en, GRB_CONTINUOUS, "dev");
		GRBVar *vars_abs = add_variables(model, en, GRB_CONTINUOUS, "abs");

		// add constraints
		add_constr_xt(model, vars_t, vars_x);
		add_constr_each_edge_dev(model, vars_t, vars_dev);
		add_constr_edge_dev_abs(model, vars_dev, vars_abs);
		//add_constr_hss(model, vars_t);//加上很差，很多模型错？

		// set objective and run model
		set_obj_min_edge_dev_abs_sum(model, vars_abs, glp_file);

		// get solution status
		set_rtn_vars_obj(model, vars_t, vars_x, rtn);
		get_overall_stats(model, lp_idx);
		set_score(model, rtn);

		// output
		if (verbose >= 2) write(cout, model, rtn);
		if (verbose == -1)
		{
			ofstream fout;
			fout.open(log_file.c_str(), ios_base::app);
			if (fout.fail()) exit(2);
			write(fout, model, rtn);
		}

		// Free environment/vars
		delete[] vars_t;
		delete[] vars_x;
		delete[] vars_dev;
		delete[] vars_abs;
		delete env;
	}
	catch (GRBException e)
	{
		cout << "Error code=" << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...)
	{
		cout << "Exception during optimization" << endl;
	}

	return rtn;// 1:have solution, 0:no solution
}

int gurobi::lp_min_edge_dev_quad_sum(string glog_file, string glp_file, int lp_idx)
{
	// minimize edge-cov-deviation-quad-sum (between the cov of an edge and the cov-sum of transcripts go through it).
	int rtn = -1; // whether got a solution
	try
	{
		GRBEnv *env = new GRBEnv(glog_file);
		GRBModel model = GRBModel(*env);
		model.set(GRB_StringAttr_ModelName, "min_edge_dev_quad_sum");
		initial_model(model);

		// add varibles
		GRBVar *vars_t = add_variables(model, tn, GRB_CONTINUOUS, "t");
		GRBVar *vars_x = add_variables(model, tn, GRB_BINARY, "x");
		GRBVar *vars_dev = add_variables(model, en, GRB_CONTINUOUS, "dev");

		// add constraints
		add_constr_xt(model, vars_t, vars_x);
		add_constr_each_edge_dev(model, vars_t, vars_dev);
		//add_constr_hss(model, vars_t);//加上很差

		// set objective and run model
		set_obj_min_edge_dev_quad_sum(model, vars_dev, glp_file);

		// get solution status
		set_rtn_vars_obj(model, vars_t, vars_x, rtn);
		get_overall_stats(model, lp_idx);
		set_score(model, rtn);

		// output
		if (verbose >= 2) write(cout, model, rtn);
		if (verbose == -1)
		{
			ofstream fout;
			fout.open(log_file.c_str(), ios_base::app);
			if (fout.fail()) exit(2);
			write(fout, model, rtn);
		}

		// Free environment/vars
		delete[] vars_t;
		delete[] vars_x;
		delete[] vars_dev;
		delete env;
	}
	catch (GRBException e)
	{
		cout << "Error code=" << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...)
	{
		cout << "Exception during optimization" << endl;
	}

	return rtn;// 1:have solution, 0:no solution
}

int gurobi::lp_hier_multiobj_dev_x(string glog_file, string glp_file, int lp_idx)
{
	// minimize max-cov-deviation and x-sum
	int rtn = -1;// whether got a solution
	try
	{
		GRBEnv *env = new GRBEnv(glog_file);
		GRBModel model = GRBModel(*env);
		model.set(GRB_StringAttr_ModelName, "lp_hier_multiobj_dev_x");
		initial_model(model);

		// add varibles
		GRBVar *vars_t = add_variables(model, tn, GRB_CONTINUOUS, "t");
		GRBVar *vars_x = add_variables(model, tn, GRB_BINARY, "x");
		GRBVar var_s = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "s");

		// add constraints
		add_constr_xt(model, vars_t, vars_x);
		add_constr_edge_dev_ub(model, vars_t, var_s);
		add_constr_hss(model, vars_t);

		// set objective and run model
		set_hier_multiobj_dev_x(model, vars_x, var_s, glp_file);

		// get solution status
		set_rtn_vars_obj(model, vars_t, vars_x, rtn);
		get_overall_stats(model, lp_idx);
		set_score(model, rtn);

		// output
		if (verbose >= 2) write(cout, model, rtn);
		if (verbose == -1)
		{
			ofstream fout;
			fout.open(log_file.c_str(), ios_base::app);
			if (fout.fail()) exit(2);
			write(fout, model, rtn);
		}

		// Free environment/vars
		delete[] vars_t;
		delete[] vars_x;
		delete env;
	}
	catch (GRBException e)
	{
		cout << "Error code=" << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...)
	{
		cout << "Exception during optimization" << endl;
	}

	return rtn;// 1:have solution, 0:no solution
}

int gurobi::lp_blend_multiobj_dev_x(string glog_file, string glp_file, int lp_idx)
{

}

int gurobi::lp_hier_multiobj_abs_x(string glog_file, string glp_file, float abs_tol, int lp_idx)
{
	// multiobj: minimize edge-cov-deviation-abs-sum and x-sum
	int rtn = -1; // whether got a solution
	try
	{
		GRBEnv *env = new GRBEnv(glog_file);
		GRBModel model = GRBModel(*env);
		model.set(GRB_StringAttr_ModelName, "lp_hier_multiobj_abs_x");
		initial_model(model);

		// add varibles
		GRBVar *vars_t = add_variables(model, tn, GRB_CONTINUOUS, "t");
		GRBVar *vars_x = add_variables(model, tn, GRB_BINARY, "x");
		GRBVar *vars_dev = add_variables(model, en, GRB_CONTINUOUS, "dev");
		GRBVar *vars_abs = add_variables(model, en, GRB_CONTINUOUS, "abs");

		// add constraints
		add_constr_xt(model, vars_t, vars_x);
		add_constr_each_edge_dev(model, vars_t, vars_dev);
		add_constr_edge_dev_abs(model, vars_dev, vars_abs);
		//add_constr_hss(model, vars_t);

		// set objective and run model
		set_hier_multiobj_abs_x(model, vars_x, vars_abs, abs_tol, glp_file);

		// get solution status
		set_rtn_vars_obj(model, vars_t, vars_x, rtn);
		get_overall_stats(model, lp_idx);
		set_score(model, rtn);

		// output
		if (verbose >= 2) write(cout, model, rtn);
		if (verbose == -1)
		{
			ofstream fout;
			fout.open(log_file.c_str(), ios_base::app);
			if (fout.fail()) exit(2);
			write(fout, model, rtn);
		}

		// Free environment/vars
		delete[] vars_t;
		delete[] vars_x;
		delete[] vars_dev;
		delete[] vars_abs;
		delete env;
	}
	catch (GRBException e)
	{
		cout << "Error code=" << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...)
	{
		cout << "Exception during optimization" << endl;
	}

	return rtn;// 1:have solution, 0:no solution
}


int gurobi::initial_model(GRBModel &model)
{
	model.set(GRB_IntParam_LogToConsole, 0);
	model.set(GRB_DoubleParam_TimeLimit, 30.0);
	model.set(GRB_IntParam_Threads, thread_count);
	model.set(GRB_DoubleParam_FeasibilityTol, 1e-9);//Primal feasibility tolerance. All constraints must be satisfied to a tolerance of FeasibilityTol.
	model.set(GRB_DoubleParam_OptimalityTol, 1e-9); //Dual feasibility tolerance. Reduced costs must all be smaller than OptimalityTol in the improving direction in order for a model to be declared optimal.
	model.set(GRB_DoubleParam_NodefileStart, 0.5); 
	if (verbose >= 3) cout << "initial_model complete.\n";
	return 0;
}

GRBVar* gurobi::add_variables(GRBModel &model, int num, char type, string name_prefix)
{
	GRBVar *vars = model.addVars(num, type);
	for (int i = 0; i < num; i++)
	{
		ostringstream ti;
		ti << name_prefix << i;
		vars[i].set(GRB_StringAttr_VarName, ti.str());
	}
	if (verbose >= 3) cout << "add_varibles complete.\n";
	return vars;
}

double* gurobi::get_coefs_x(double exponent)
{
	double *coefs = new double[tn];
	for (int i = 0; i < tn; i++)
	{
		int a = 0;
		for (int j = 0; j < hsn; j++) a += hss[j][i];
		if (a <= 0) a = 1;
		coefs[i] = pow(a, -3);
	}
	return coefs;
}

int gurobi::add_constr_xt(GRBModel &model, GRBVar *vars_t, GRBVar *vars_x)
{
	for (int i = 0; i < tn; i++)
	{
		model.addConstr(vars_x[i] <= vars_t[i]*ALPHA);
		model.addConstr(vars_t[i]*coef <= vars_x[i]*BETA);
	}
	if (verbose >= 3) cout << "add_constr_xt complete." << endl;
	return 0;
}

int gurobi::add_constr_edge_dev_ub(GRBModel &model, GRBVar *vars_t, GRBVar &var_s)
{
	for (int i = 0; i < en; i++)
	{
		GRBLinExpr lhs1 = 0, lhs2 = 0;
		for (int j = 0; j < tn; j++)
		{
			lhs1 += vars_t[j] * B[j][i];
			lhs2 += vars_t[j] * B[j][i];
		}
		lhs1 -= var_s;
		lhs2 += var_s;
		model.addConstr(lhs1 <= ew[i]);
		model.addConstr(lhs2 >= ew[i]);
	}
	if (verbose >= 3) cout << "add_constr_edge_dev_ub complete." << endl;
	return 0;
}

int gurobi::add_constr_edge_dev_ub(GRBModel &model, GRBVar *vars_t, double max_dev)
{
	for (int i = 0; i < en; i++)
	{
		GRBLinExpr lhs1 = 0, lhs2 = 0;
		for (int j = 0; j < tn; j++)
		{
			lhs1 += vars_t[j] * B[j][i];
			lhs2 += vars_t[j] * B[j][i];
		}
		lhs1 -= max_dev;
		lhs2 += max_dev;
		model.addConstr(lhs1 <= ew[i]);
		model.addConstr(lhs2 >= ew[i]);
	}
	if (verbose >= 3) cout << "add_constr_edge_dev_ub complete." << endl;
	return 0;
}

int gurobi::add_constr_each_edge_dev(GRBModel &model, GRBVar *vars_t, GRBVar *vars_dev)
{
	for (int i = 0; i < en; i++)
	{
		GRBLinExpr lhs = 0;
		lhs += ew[i];
		for (int j = 0; j < tn; j++) lhs -= vars_t[j] * B[j][i];
		model.addConstr(vars_dev[i] == lhs);
	}
	if (verbose >= 3) cout << "add_constr_each_edge_dev complete." << endl;
	return 0;
}

int gurobi::add_constr_edge_dev_abs(GRBModel &model, GRBVar *vars_dev, GRBVar *vars_abs)
{
	for (int i = 0; i < en; i++) model.addGenConstrAbs(vars_abs[i], vars_dev[i]);
	if (verbose >= 3) cout << "add_constr_edge_dev_abs complete." << endl;
	return 0;
}

int gurobi::add_constr_hss(GRBModel &model, GRBVar *vars_t)
{
	for (int i = 0; i < hss.size(); i++)
	{
		GRBLinExpr lhs = 0;
		for (int j = 0; j < tn; j++)
		{
			lhs += vars_t[j] * hss[i][j];
		}
		model.addConstr(lhs >= ucnts[i]);
	}
	if (verbose >= 3) cout << "add_constr_hss complete." << endl;
	return 0;
}


int gurobi::set_obj_min_edge_diff(GRBModel &model, GRBVar &var_s, string glp_file)
{
	GRBLinExpr obj = 0;
	obj += var_s;
	model.setObjective(obj, GRB_MINIMIZE);
	model.optimize();
	model.write(glp_file);
	if (verbose >= 3) cout << "set_obj_min_edge_diff complete." << endl;
	return 0;
}

int gurobi::set_obj_min_x_sum(GRBModel &model, GRBVar *vars_x, string glp_file)
{
	GRBLinExpr obj = 0;
	for (int i = 0; i < tn; i++)
	{
		int a = 0;
		for (int j = 0; j < hsn; j++) a += hss[j][i];
		if (a <= 0) a = 1;
		double xcoef = 1 / (double)(a*a*a);
		obj += xcoef*vars_x[i];
	}
	model.setObjective(obj, GRB_MINIMIZE);
	model.optimize();
	model.write(glp_file);
	if (verbose >= 3) cout << "set_obj_min_x_sum complete." << endl;
	return 0;
}

int gurobi::set_obj_min_edge_dev_abs_sum(GRBModel &model, GRBVar *vars_abs, string glp_file)
{
	GRBLinExpr lhs = 0;
	for (int i = 0; i < en; i++) lhs += vars_abs[i];
	model.setObjective(lhs, GRB_MINIMIZE);
	model.optimize();
	model.write(glp_file);
	if (verbose >= 3) cout << "set_obj_min_edge_dev_abs_sum complete." << endl;
	return 0;
}

int gurobi::set_obj_min_edge_dev_quad_sum(GRBModel &model, GRBVar *vars_dev, string glp_file)
{
	GRBQuadExpr qexpr;
	double *coefs = new double[en];
	for (int i = 0; i < en; i++) coefs[i] = 1;
	qexpr.addTerms(coefs, vars_dev, vars_dev, en);
	model.setObjective(qexpr, GRB_MINIMIZE);
	model.optimize();
	model.write(glp_file);
	if (verbose >= 3) cout << "set_obj_min_edge_dev_quad_sum complete." << endl;
	delete[] coefs;
	return 0;
}

int gurobi::set_hier_multiobj_dev_x(GRBModel &model, GRBVar *vars_x, GRBVar &var_s, string glp_file)
{
	model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

	GRBLinExpr obj0 = var_s;
	model.setObjectiveN(obj0, 0, 1, 1, 0, 0, "min_max_dev");

	GRBLinExpr obj1;
	double *coefs = get_coefs_x(x_coef_exp);
	obj1.addTerms(coefs, vars_x, tn);
	model.setObjectiveN(obj1, 1, 0, 1, 0, 0, "min_x_sum");
	model.optimize();
	model.write(glp_file);
	if (verbose >= 3) cout << "set_hier_multiobj_dev_x complete." << endl;
	delete[] coefs;
	return 0;
}

int gurobi::set_hier_multiobj_abs_x(GRBModel &model, GRBVar *vars_x, GRBVar *vars_abs, float reltol, string glp_file)
{
	model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

	GRBLinExpr obj0;
	double *coefs0 = new double[en];
	for (int i = 0; i < en; i++) coefs0[i] = 1;
	obj0.addTerms(coefs0, vars_abs, en);
	model.setObjectiveN(obj0, 0, 1, 1, 0, reltol, "min_abs_sum");

	GRBLinExpr obj1;
	double *coefs1 = get_coefs_x(x_coef_exp);
	obj1.addTerms(coefs1, vars_x, tn);
	model.setObjectiveN(obj1, 1, 0, 1, 0, 0, "min_x_sum");

	model.optimize();
	model.write(glp_file);
	if (verbose >= 3) cout << "set_hier_multiobj_abs_x complete." << endl;
	delete[] coefs0;
	delete[] coefs1;
	return 0;
}

// set rtn: 2:opt, 1:tl -1:no-solution
int gurobi::set_rtn_vars_obj(GRBModel &model, GRBVar *vars_t, GRBVar *vars_x, int &rtn)
{
	int status = model.get(GRB_IntAttr_Status);
	if (status == GRB_INF_OR_UNBD) rtn = -1;
	else if (status == GRB_INFEASIBLE) rtn = -1;
	else if (status == GRB_UNBOUNDED) rtn = -1;
	else if (status == GRB_OPTIMAL)
	{
		rtn = 2;
		ts.clear();
		xs.clear();
		for (int i = 0; i < tn; i++)
		{
			ts.push_back(vars_t[i].get(GRB_DoubleAttr_X));
			xs.push_back(vars_x[i].get(GRB_DoubleAttr_X));
		}
		if (verbose >= 4)
		{
			double obj_value = model.get(GRB_DoubleAttr_ObjVal);
			cout << "ObjVal=" << obj_value << endl;
		}
		int objcount = model.get(GRB_IntAttr_NumObj);
		if (objcount == 1) obj_values.push_back(model.get(GRB_DoubleAttr_ObjVal));
		else if (objcount > 1)
		{
			for (int i = 0; i < objcount; i++)
			{
				model.set(GRB_IntParam_ObjNumber, i);
				obj_values.push_back(model.get(GRB_DoubleAttr_ObjNVal));
			}
		}
		if (verbose >= 4)
		{
			cout << "obj_values: ";
			writev(cout, obj_values);
			cout << endl;
		}
	}
	else if (status == GRB_TIME_LIMIT)
	{
		int solcount = model.get(GRB_IntAttr_SolCount);
		if (solcount <= 0) rtn = -1;
		else
		{
			rtn = 1;
			ts.clear();
			xs.clear();
			model.set(GRB_IntParam_SolutionNumber, 0);
			for (int i = 0; i < tn; i++)
			{
				ts.push_back(vars_t[i].get(GRB_DoubleAttr_Xn));
				xs.push_back(vars_x[i].get(GRB_DoubleAttr_Xn));
			}
			if (verbose >= 4)
			{
				double obj_value = model.get(GRB_DoubleAttr_PoolObjVal);
				cout << "PoolObjVal=" << obj_value << endl;
			}
			int objcount = model.get(GRB_IntAttr_NumObj);
			if (objcount == 1) obj_values.push_back(model.get(GRB_DoubleAttr_ObjVal));
			else if (objcount > 1)
			{
				for (int i = 0; i < objcount; i++)
				{
					model.set(GRB_IntParam_ObjNumber, i);
					obj_values.push_back(model.get(GRB_DoubleAttr_ObjNVal));
				}
			}
			if (verbose >= 4)
			{
				cout << "obj_values: ";
				writev(cout, obj_values);
				cout << endl;
			}
		}
	}
	else rtn = -1;
	return 0;
}

int gurobi::get_overall_stats(GRBModel &model, int lp_idx)
{
	int status = model.get(GRB_IntAttr_Status);
	lp_idx == 1 ? lp1_stats[0] += 1 : lp2_stats[0] += 1;
	int solcount = model.get(GRB_IntAttr_SolCount);
	switch (status)
	{
	case GRB_INF_OR_UNBD:
		if (lp_idx == 1)
		{
			lp1merr_gids += (gid + " ");
			lp1_stats[1] += 1;
		}
		else
		{
			lp2merr_gids += (gid + " ");
			lp2_stats[1] += 1;
		}
		break;
	case GRB_INFEASIBLE:
		if (lp_idx == 1)
		{
			lp1merr_gids += (gid + " ");
			lp1_stats[1] += 1;
		}
		else
		{
			lp2merr_gids += (gid + " ");
			lp2_stats[1] += 1;
		}
		break;
	case GRB_UNBOUNDED:
		if (lp_idx == 1)
		{
			lp1merr_gids += (gid + " ");
			lp1_stats[1] += 1;
		}
		else
		{
			lp2merr_gids += (gid + " ");
			lp2_stats[1] += 1;
		}
		break;
	case GRB_OPTIMAL:
		lp_idx == 1 ? lp1_stats[2] += 1 : lp2_stats[2] += 1;
		break;
	case GRB_TIME_LIMIT:
		if (solcount <= 0)
		{
			if (lp_idx == 1)
			{
				lp1nsol_gids += (gid + " ");
				lp1_stats[4] += 1;
			}
			else
			{
				lp2nsol_gids += (gid + " ");
				lp2_stats[4] += 1;
			}
		}
		else lp_idx == 1 ? lp1_stats[3] += 1 : lp2_stats[3] += 1;
		break;
	default:
		lp_idx == 1 ? lp1_stats[5] += 1 : lp2_stats[5] += 1;
		break;
	}
	return 0;
}

int gurobi::set_score(GRBModel &model, int rtn)
{
	score = rtn;
	if (score == -1) return 0;
	double max_ew;
	switch(lp_option)
	{
	case 4:
		//lpHierDevX
		assert(obj_values.size() == 2);
		max_ew = *max_element(ew.begin(), ew.end());
		if (obj_values[0] > max_ew) score = 0;
		break;
	case 6:
		//lpHierAbsX
		assert(obj_values.size() == 2);
		break;
	default:
		break;
	}
	return 0;
}

int gurobi::write(ostream &out, GRBModel &model, int rtn)
{
	write_lp_status(out, model);
	if (rtn >= 0) write_t_x_obj(out, model);
	return 0;
}

int gurobi::write_lp_status(ostream &out, GRBModel &model)
{
	string name = model.get(GRB_StringAttr_ModelName);
	out << "model=" << name << ", ";
	int status = model.get(GRB_IntAttr_Status);
	if (status == GRB_INF_OR_UNBD)
	{
		out << "status = GRB_INF_OR_UNBD = " << status << "\n";
		return -1;
	}
	else if (status == GRB_INFEASIBLE)
	{
		out << "status = GRB_INFEASIBLE = " << status << ", get no solution.\n";
		// do IIS
		out << "Computing IIS, the following constraint(s) cannot be satisfied:\n";
		model.computeIIS();
		GRBConstr* c = 0;
		c = model.getConstrs();
		for (int i = 0; i < model.get(GRB_IntAttr_NumConstrs); ++i)
		{
			if (c[i].get(GRB_IntAttr_IISConstr) == 1)
			{
				out << c[i].get(GRB_StringAttr_ConstrName) << " ";
			}
		}
		out << "\n";
		return -1;
	}
	else if (status == GRB_UNBOUNDED)
	{
		out << "status = GRB_UNBOUNDED = " << status << ", get no solution.\n";
		return -1;
	}
	else if (status == GRB_OPTIMAL)
	{
		out << "status = GRB_OPTIMAL = " << status << ", get optimal solution:\n";
		return 0;
	}
	else if (status == GRB_TIME_LIMIT)
	{
		out << "status = GRB_TIME_LIMIT = " << status << ", ";
		int solcount = model.get(GRB_IntAttr_SolCount);
		if (solcount <= 0)
		{
			out << "get no solution.\n";
			return -1;
		}
		else
		{
			out << "get feasible solution:\n";
			return 0;
		}
	}
	else
	{
		out << "Into weird status: " << status << ", get no solution.\n";
		return -1;
	}
}

int gurobi::write_t_x_obj(ostream &out, GRBModel &model)
{
	assert(ts.size() == tn);
	assert(xs.size() == tn);
	int objcount = model.get(GRB_IntAttr_NumObj);
	if (objcount == 1) out << "obj=" << model.get(GRB_DoubleAttr_ObjVal) << "\n";
	else if (objcount > 1)
	{
		for (int i = 0; i < objcount; i++)
		{
			model.set(GRB_IntParam_ObjNumber, i);
			out << "obj" << i << "=" << model.get(GRB_StringAttr_ObjNName) << "=" << model.get(GRB_DoubleAttr_ObjNVal) << '\n';
		}
	}
	out << "final solution:\nobj=";
	writev(out, obj_values);
	out << "\n";
	for (int i = 0; i < ts.size(); i++)
	{
		out << i << "={" << xs[i] << "," << ts[i] << "} ";
	}
	out << "\n";
	return 0;
}
