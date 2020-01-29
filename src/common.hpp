#ifndef COMMON_H_
#define COMMON_H_

#include <boost/program_options.hpp>
#include <iostream>
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

using namespace std;
namespace po = boost::program_options;

#define VERBOSE 2                   // printing level  (=10 only incumbent, =20 little output, =50-60 good, =70 verbose, >=100 cplex log)

//hard-wired parameters
#define EPSILON_MED 1e-5
#define EPSILON_SMALL 1e-9
#define TOLERANCE 1e-1

#define LAZY_CONST_MODEL_TYPE 1
#define NEW_ROWS_MODEL_TYPE 4
#define CALLBACK_MODEL_TYPE 3
#define LOOP_TASK_MODEL_TYPE 2

// data structures
typedef struct
{
    // input data
    int n_turbines;
    int n_cables;
    int subs_cables;					// max number of input cables acceptable by the substation
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
    bool print_model;
    bool rins;
    string solver;
    int model_type;
    bool slack_substation;				// add slack variable to in-deg substation constraint (in order to find quicly an incumbent solution)
    bool slack_flow;					// add slack variables to flow constraints
    int randomseed;
    int num_threads;
    double time_limit;					// overall time limit, in seconds
    double elapsed_time;
    string turbine_input_file;			// turbine input file (name)
    string cable_input_file;			// cable input file (name)
    double iter_time;					// iteration time
    double hard_fix_prob;				// hard fix probability
    bool use_cross_table;
    bool multi; 						// use multi-threads optimisation when possible
    int ncores;

    // global data
    struct timespec *time_start;
    struct timespec *time_end;
    char time_limit_expired;
    double tbest;						// time for the best sol. available
    double *xstar;						// best solution available
    double best_lb;						// best lower bound available
    double obj_val;

    // callback data structures
    double **xstar_vector;				// poiters to xstar vectors for thread safe mode
    int **indeces_vector;
    double **values_vector;

    // model
    int fstart;
    int xstart;
    int ystart;
    int sstart;
    int lstart;
    int ncols;							// number of columns
    int indeg_constraint;				// row that contains the index of indeg constraint
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

inline int f_pos(int i, int j, instance *inst) { return inst->fstart + i * inst->n_turbines + j; }
inline int x_pos(int i, int j, int k, instance *inst) { return inst->xstart + i * inst->n_turbines * inst->n_cables + j * inst->n_cables + k; }
inline int y_pos(int i, int j, instance *inst) { return inst->ystart + i * inst->n_turbines + j; }

int comm_free_instance(instance *inst);
int comm_read_input(instance *inst);

int comm_is_crossing(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);
int comm_no_cross_separation(instance *inst, double *xstar, int i, int j, int k, int h, int thread, int table);

char*	comm_generate_cross_table(instance *inst);
void	comm_set_cross_table(instance *inst, char *cross_table, char value, int i, int j, int k, int h);
int		comm_get_cross_table(instance *inst, char *cross_table, int i, int j, int k, int h);

void	comm_print_error(const char *err);
int		comm_is_time_limit_expired(instance *inst);
double	comm_get_time_elapsed(instance *inst);

#endif //COMMON_H_
