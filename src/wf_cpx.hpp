
#ifndef WF_CPX_H_
#define WF_CPX_H_

#include "common.hpp"
#include <cplex.h>

#define TICKS_PER_SECOND 1000.0 // cplex's ticks on Intel Core i7 quadcore @2.3GHZ

int WFopt(instance *inst);
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp);
void build_model_cross_constraints(instance *inst, CPXENVptr env, CPXLPptr lp);
void loop_solver(instance *inst, CPXENVptr env, CPXLPptr lp);
void callback_solver(instance *inst, CPXENVptr env, CPXLPptr lp);
void hardfix_solver(instance *inst, CPXENVptr env, CPXLPptr lp);

int compute_no_cross_cuts(instance *inst, CPXCENVptr env, CPXLPptr lp, int thread);
int CPXPUBLIC callback_cross_multi(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p);
int CPXPUBLIC callback_cross_single(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p);

int hard_fix(instance *inst, CPXENVptr env, CPXLPptr lp);
int unfix_variables(instance *inst, CPXENVptr env, CPXLPptr lp);

int mip_solved_to_optimality(CPXENVptr env, CPXLPptr lp);

#endif //WF_CPX_H_
