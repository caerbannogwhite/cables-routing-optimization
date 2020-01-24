#ifndef AUX_H_
#define AUX_H_

#include "wf.h"

int is_crossing(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);
int no_cross_separation(instance *inst, double *xstar, int i, int j, int k, int h, int thread, int table);
int compute_no_cross_cuts(instance *inst, CPXCENVptr env, CPXLPptr lp, int thread);
int CPXPUBLIC callback_cross_multi(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p);
int CPXPUBLIC callback_cross_single(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p);

char *generate_cross_table(instance *inst);
void set_cross_table(instance *inst, char *cross_table, char value, int i, int j, int k, int h);
int get_cross_table(instance *inst, char *cross_table, int i, int j, int k, int h);

int hard_fix(instance *inst, CPXENVptr env, CPXLPptr lp);
int unfix_variables(instance *inst, CPXENVptr env, CPXLPptr lp);

#endif //AUX_H_
