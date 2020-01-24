#ifndef WF_H_
#define WF_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <cplex.h>
#include <pthread.h>
#include <time.h>

#define VERBOSE				  	50		// printing level  (=10 only incumbent, =20 little output, =50-60 good, =70 verbose, >=100 cplex log)

//hard-wired parameters
#define EPSILON_MED				1e-5
#define EPSILON_SMALL			1e-9		// 1e-9		// very small numerical tolerance
#define TOLERANCE 				1e-1
#define TICKS_PER_SECOND		1000.0  	// cplex's ticks on Intel Core i7 quadcore @2.3GHZ

#define LAZY_CONST_MODEL_TYPE   1
#define NEW_ROWS_MODEL_TYPE     4
#define CALLBACK_MODEL_TYPE     3
#define LOOP_TASK_MODEL_TYPE    2

#define DBL_MAX 10e20 /* max value */

// data structures
typedef struct {

    // input data
    int n_turbines;
    int n_cables;
    int subs_cables;						// max number of input cables for the substation
    int subs_index;
    double *x_turb_coords;
    double *y_turb_coords;
    double *turb_powers;
    double turb_tot_power;
    double x_subs_coord;
    double y_subs_coord;
    double *cable_powers;
    double *cable_costs;

    // parameters
    int print_model;
    int rins;
    double polishafter;
    char solver[1000];
    int model_type;
    int slack_substation;					// add slack variable to in-deg substation constraint (in order to find quicly an incumbent solution)
    int slack_flow;                         // add slack variables to flow constraints
    int randomseed;
    int num_threads;
    double time_limit;						// overall time limit, in sec.s
    double elapsed_time;
    char turbine_input_file[1000];	  		// turbine input file (name)
    char cable_input_file[1000];			// cable input file (name)
    char node_file[1000];		  			// cplex node file
    char params_file[1000];
    double iter_time;                       // iteration time
    double hard_fix_prob;                   // hard fix probability
    int use_cross_table;
    int multi;                              // use multi-threads if possible
    int ncores;

    // global data
    struct timespec *time_start;
    struct timespec *time_end;
    char time_limit_expired;
    double tbest;							// time for the best sol. available
    double *xstar;							// best sol. available
    double best_lb;							// best lower bound available
    double obj_val;

    // callback data structures
    double **xstar_vector;					// poiters to xstar vectors for thread safe mode
    int **indeces_vector;
    double **values_vector;

    // model
    int fstart;
    int xstart;
    int ystart;
    int sstart;
    int lstart;
    int ncols;								// number of colons
    int indeg_constraint; 					// row that contains the index of indeg constraint
    int flow_constraints_start;

    // heur data structures
    int *preds;
    int *preds_copy;
    int *flags;

    double *solution;
    double *streams;

    char *cross_table;

    // debug
    int debug;
    int counter;
    int tot_time;
} instance;

//inline
inline int imax(int i1, int i2) { return i1 > i2 ? i1 : i2; }
inline int imin(int i1, int i2) { return i1 < i2 ? i1 : i2; }
inline double dmax(double d1, double d2) { return d1 > d2 ? d1 : d2; }
inline double dmin(double d1, double d2) { return d1 < d2 ? d1 : d2; }

inline int fpos(int i, int j, instance *inst) { return inst->fstart + i * inst->n_turbines + j; }
inline int xpos(int i, int j, int k, instance *inst) { return inst->xstart + i * inst->n_turbines * inst->n_cables + j * inst->n_cables + k; }
inline int ypos(int i, int j, instance *inst) { return inst->ystart + i * inst->n_turbines + j; }

int WFopt(instance *inst);
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp);
void build_model_cross_constraints(instance *inst, CPXENVptr env, CPXLPptr lp);
void loop_solver(instance *inst, CPXENVptr env, CPXLPptr lp);
void callback_solver(instance *inst, CPXENVptr env, CPXLPptr lp);
void hardfix_solver(instance *inst, CPXENVptr env, CPXLPptr lp);

void print_error(const char *err);
char *get_log();
int is_time_limit_expired(instance *inst);
double get_time_elapsed(instance *inst);
int mip_solved_to_optimality(CPXENVptr env, CPXLPptr lp);

#endif //WF_H_
