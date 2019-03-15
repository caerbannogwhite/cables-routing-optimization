#include "wf.h"
#include "aux.h"
#include "heur.h"

int WFopt(instance *inst);
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp);
void build_model_cross_constraints(instance *inst, CPXENVptr env, CPXLPptr lp);
void loop_solver(instance *inst, CPXENVptr env, CPXLPptr lp);
void callback_solver(instance *inst, CPXENVptr env, CPXLPptr lp);
void hardfix_solver(instance *inst, CPXENVptr env, CPXLPptr lp);

int WFopt(instance *inst)
{
    int error, i, j, k, ncols;
    FILE *fout;

    // BUILD MODEL AND SOLVE
    inst->ystart = -1;
    inst->fstart = -1;
    inst->xstart = -1;
    inst->sstart = -1;

    printf("%s Start building model.\n", get_log());
    if (!(strcmp(inst->solver, "VNS") == 0 || strcmp(inst->solver, "simulated_annealing") == 0))
    {
        // OPEN CPLEX
        CPXENVptr env = CPXopenCPLEX(&error);
        CPXLPptr lp = CPXcreateprob(env, &error, "WF");

        // SETTING PARAMETERS
        CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 5);               // cplex output
        CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);              // cplex logs

        CPXsetintparam(env, CPX_PARAM_RINSHEUR, inst->rins);        // rins heuristic
        CPXsetdblparam(env, CPX_PARAM_TILIM, inst->time_limit);     // cplex mipopt time limit

        // precision
        CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0);
        //CPXsetdblparam(env, CPX_PARAM_EPGAP, 1e-5);
        //CPXsetdblparam(env, CPX_PARAM_EPAGAP, 1e-2);
        CPXsetdblparam(env, CPX_PARAM_EPRHS, 1e-9);

        // initial LP algorithm
        //CPXsetintparam(env, CPX_PARAM_LPMETHOD, 4);   // barrier at the root

        // barrier vs simplex
        //CPXsetintparam(env, CPX_PARAM_STARTALG, 6);   // concurrent, cplex bug??
        
        // polish
        // TODO: sistemare parametri
        CPXsetdblparam(env, CPX_PARAM_POLISHAFTERDETTIME, inst->polishafter*TICKS_PER_SECOND);
        //CPXsetintparam(env, CPX_PARAM_POLISHAFTERNODE, 0);
        //CPXsetdblparam(env, CPX_PARAM_DETTILIM, xxx*TICKS_PER_SECOND+0.0);

        // random seed
        if (inst->randomseed != 0) CPXsetintparam(env, CPX_PARAM_RANDOMSEED, abs(inst->randomseed));

        if (strcmp(inst->solver, "none") == 0)
        {
            build_model(inst, env, lp);
            if (inst->print_model) CPXwriteprob(env, lp, "model.lp", NULL);
        
            CPXmipopt(env, lp);

            inst->xstar = (double *) calloc(sizeof(double), inst->ncols);
            CPXgetx(env, lp, inst->xstar, 0, inst->ncols-1);
        }
        else if (strcmp(inst->solver, "lazy") == 0)
        {
            inst->model_type = LAZY_CONST_MODEL_TYPE;
            build_model_cross_constraints(inst, env, lp);
            if (inst->print_model) CPXwriteprob(env, lp, "model.lp", NULL);
        
            CPXmipopt(env, lp);

            inst->xstar = (double *) calloc(inst->ncols, sizeof(double));
            CPXgetx(env, lp, inst->xstar, 0, inst->ncols-1);
        }
        else if (strcmp(inst->solver, "newrows") == 0)
        {
            inst->model_type = NEW_ROWS_MODEL_TYPE;
            build_model_cross_constraints(inst, env, lp);
            if (inst->print_model) CPXwriteprob(env, lp, "model.lp", NULL);
        
            CPXmipopt(env, lp);

            inst->xstar = (double *) calloc(inst->ncols, sizeof(double));
            CPXgetx(env, lp, inst->xstar, 0, inst->ncols-1);
        }
        else if (strcmp(inst->solver, "loop") == 0)
        {
            // build the default model
            build_model(inst, env, lp);
            inst->xstar = (double *) calloc(inst->ncols, sizeof(double));
            loop_solver(inst, env, lp);
        }
        else if (strcmp(inst->solver, "callback") == 0)
        {
            callback_solver(inst, env, lp);
            inst->xstar = (double *) calloc(inst->ncols, sizeof(double));
            if (CPXgetx(env, lp, inst->xstar, 0, inst->ncols-1)) printf("%s Error: CPXgetx no solution. \n", get_log());;
        }
        else if (strcmp(inst->solver, "hardfix") == 0)
        {
            inst->model_type = LAZY_CONST_MODEL_TYPE;
            build_model_cross_constraints(inst, env, lp);
            inst->xstar = (double *) calloc(inst->ncols, sizeof(double));
            hardfix_solver(inst, env, lp);
        }
        else
        {
            print_error("Solver unknown.");
        }

        if (inst->print_model) CPXwriteprob(env, lp, "model.lp", NULL);

        //if (strcmp(inst->solver, "hardfix") == 0)
        //{
        //    CPXgetbestobjval(env, lp, &inst->best_lb);
        //    printf("\n%s Best obj-val=%lf\n", get_log(), inst->best_lb);
        //    CPXgetobjval(env, lp, &inst->obj_val);
        //    printf("%s Best integer=%lf\n", get_log(), inst->obj_val);
        //}

        // free pools and close cplex model
        printf("%s Closing CPLEX.\n", get_log());
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
    }
    else if (strcmp(inst->solver, "VNS") == 0)
    {
        heur_VNS_launcher(inst);
    }
    else if (strcmp(inst->solver, "simulated_annealing") == 0)
    {
        heur_simulated_annealing_launcher(inst);
    }
    else
    {
        print_error("Solver unknown.");
    }

    // PRINT SOLUTION
    fout = fopen("solution", "w");
    for (i = 0; i < inst->n_turbines; i++)
    {
        for (j = 0; j < inst->n_turbines; j++)
            if (inst->xstar[fpos(i, j, inst)] > EPSILON_SMALL) fprintf(fout, "f %d %d %d %lf\n", i, j, -1, inst->xstar[fpos(i, j, inst)]);
    }

    for (i = 0; i < inst->n_turbines; i++)
    {
        for (j = 0; j < inst->n_turbines; j++)
        {
            for (k = 0; k < inst->n_cables; k++)
                if (inst->xstar[xpos(i, j, k, inst)] > EPSILON_SMALL) fprintf(fout, "x %d %d %d %lf\n", i, j, k, inst->xstar[xpos(i, j, k, inst)]);
        }
    }

    for (i = 0; i < inst->n_turbines; i++)
    {
        for (j = 0; j < inst->n_turbines; j++)
            if (inst->xstar[ypos(i, j, inst)] > EPSILON_SMALL) fprintf(fout, "y %d %d %d %lf\n", i, j, -1, inst->xstar[ypos(i, j, inst)]);
    } 
    fclose(fout);

    return 0;
}

/*************************************************************************************************/
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp)
/*************************************************************************************************/
{
    char binary = 'B';
    char continuous = 'C';
    char integer = 'I';
    char equal = 'E';
    char greater_equal = 'G';
    char less_equal = 'L';
    char **cname;
    int h, i, j, k, lastrow;
    double obj;
    double obj_slack = 1e12; // obj function value of slack variables
    double zero = 0.0;
    double one = 1.0;
    double up = 0.0;
    double cable_max_value = inst->cable_powers[inst->n_cables - 1];

    cname = (char **) calloc(1, sizeof(char *));		// (char **) required by cplex
    cname[0] = (char *) calloc(100, sizeof(char));


    // COLUMNS
    // add general integer var f(i,j): value of current flow in (i,j)
    if (inst->fstart != -1) print_error("Error in build_model0(): var f cannot be redefined.");
    inst->fstart = CPXgetnumcols(env,lp);
    printf("%s Adding f variables to the model at column %d.\n", get_log(), inst->fstart);
    for (i = 0; i < inst->n_turbines; i++)
    {
        for (j = 0; j < inst->n_turbines; j++)
        {
            // skip loop and double edges
            up = (i == j) ? 0.0 : CPX_INFBOUND;

            //sprintf(cname[0], "f(%d,%d)", i+1, j+1);
            sprintf(cname[0], "f(%d,%d)", i, j);

            // env: CPLEX environment
            // lp: CPLEX problem object
            // ccnt: number of new variables being added (always 1 here)
            // obj: objective function coefficient
            // lb: lower bound
            // ub: upper bound, NULL means CPX_INFBOUND
            // xctype: general integer
            // colname: name of the variable
            if (CPXnewcols(env, lp, 1, NULL, &zero, &up, &continuous, cname)) print_error(" wrong CPXnewcols on f vars.");
            if (CPXgetnumcols(env,lp)-1 != fpos(i, j, inst)) print_error(" wrong position for f vars.");
        }
    }

    // add binary var x(i,j,k): in connecting (i,j) a k-value cable is used
    if (inst->xstart != -1) print_error("Error in build_model0(): var x cannot be redefined.");
    inst->xstart = CPXgetnumcols(env,lp);
    printf("%s Adding x variables to the model at column %d.\n", get_log(), inst->xstart);
    for (i = 0; i < inst->n_turbines; i++)
    {
        for (j = 0; j < inst->n_turbines; j++)
        {
            // skip loop and double edges
            up = (i == j) ? 0.0 : 1.0;

            for (k = 0; k < inst->n_cables; k++)
            {
                //sprintf(cname[0], "x(%d,%d,%d)", i+1, j+1, k+1);
                sprintf(cname[0], "x(%d,%d,%d)", i, j, k);
                obj = inst->cable_costs[k] * hypot(inst->x_turb_coords[i] - inst->x_turb_coords[j], inst->y_turb_coords[i] - inst->y_turb_coords[j]);

                // obj: objective function coefficient, cost(x) * dist(x)
                if (CPXnewcols(env, lp, 1, &obj, &zero, &up, &binary, cname)) print_error(" wrong CPXnewcols on x vars.");
                if (CPXgetnumcols(env,lp)-1 != xpos(i, j, k, inst)) print_error(" wrong position for x vars.");
            }
        }
    }

    // add binary var y(i,j): connection in (i,j)
    if (inst->ystart != -1) print_error("Error in build_model0(): var y cannot be redefined.");
    inst->ystart = CPXgetnumcols(env,lp);
    printf("%s Adding y variables to the model at column %d.\n", get_log(), inst->ystart);
    for (i = 0; i < inst->n_turbines; i++)
    {
        for (j = 0; j < inst->n_turbines; j++)
        {
            // skip loop and double edges
            up = (i == j) ? 0.0 : 1.0;

            //sprintf(cname[0], "y(%d,%d)", i+1, j+1);
            sprintf(cname[0], "y(%d,%d)", i, j);

            if (CPXnewcols(env, lp, 1, NULL, &zero, &up, &binary, cname)) print_error(" wrong CPXnewcols on y vars.");
            if (CPXgetnumcols(env,lp)-1 != ypos(i, j, inst)) print_error(" wrong position for y vars.");
        }
    }

    // substation in-deg slack variable
    if (inst->slack_substation)
    {
        inst->sstart = CPXgetnumcols(env,lp);
        printf("%s Adding s variable at column %d.\n", get_log(), inst->sstart);
        sprintf(cname[0], "s(0)");
        if (CPXnewcols(env, lp, 1, &obj_slack, &zero, NULL, &integer, cname)) print_error("Wrong CPXnewcols() on slack variables.");
    }

    // flow slack variables
    if (inst->slack_flow)
    {
        inst->lstart = CPXgetnumcols(env,lp);
        printf("%s Adding l variables at column %d.\n", get_log(), inst->lstart);

        for (i = 0; i < inst->n_turbines; ++i)
        {
            sprintf(cname[0], "l(%d)", i);
            if (CPXnewcols(env, lp, 1, &obj_slack, &zero, &cable_max_value, &continuous, cname)) print_error("Wrong CPXnewcols() on slack variables.");
        }
    }


    // ROWS
    // turbine (and substation) out-degree constraints
    printf("%s Adding turbine out-degree constraints at row %d.\n", get_log(), CPXgetnumrows(env,lp));
    for (h = 0; h < inst->n_turbines; ++h)
    {
        lastrow = CPXgetnumrows(env,lp);
        //sprintf(cname[0], "outdeg(%d)", h+1);
        sprintf(cname[0], "outdeg(%d)", h);

        // env: CPLEX environment
        // lp: CPLEX problem object
        // rcnt: number of new rows being added (always 1 here)
        // rhs: right hand side
        // sense: 'E' for equations
        // rngval: range values for the new constraints
        // rowname: name of the constraint

        // looking for the substation
        if (inst->turb_powers[h] < -0.5)
        {
            if (CPXnewrows(env, lp, 1, &zero, &equal, NULL, cname)) print_error("Wrong CPXnewrows() on turbine out-degree.");
        }
        else
        {
            if (CPXnewrows(env, lp, 1, &one, &equal, NULL, cname)) print_error("Wrong CPXnewrows() on turbine out-degree.");
        }
        for (j = 0; j < inst->n_turbines; j++)
            if (CPXchgcoef(env, lp, lastrow, ypos(h, j, inst), 1.0)) print_error("Wrong CPXchgcoef() on turbine out-degree.");
    }

    // substation in-degree constraint
    printf("%s Adding substation in-degree constraint at row %d.\n", get_log(), CPXgetnumrows(env,lp));
    for (h = 0; h < inst->n_turbines; ++h)
    {
        if (inst->turb_powers[h] < -0.5)
        {
            lastrow = CPXgetnumrows(env,lp);
            inst->indeg_constraint = lastrow;
            //sprintf(cname[0], "subs_indeg(%d)", h+1);
            sprintf(cname[0], "subs_indeg(%d)", h);
            double rhs = (double) inst->subs_cables;

            if (CPXnewrows(env, lp, 1, &rhs, &less_equal, NULL, cname)) print_error("Wrong CPXnewrows() on substation in-degree.");
            for (i = 0; i < inst->n_turbines; i++)
                if (CPXchgcoef(env, lp, lastrow, ypos(i, h, inst), 1.0)) print_error("Wrong CPXchgcoef() on substation in-degree.");
            break;
        }
    }

    // turbine flow constraints
    printf("%s Adding turbine flow constraints at row %d.\n", get_log(), CPXgetnumrows(env,lp));
    inst->flow_constraints_start = CPXgetnumrows(env, lp);
    for (h = 0; h < inst->n_turbines; ++h)
    {

        lastrow = CPXgetnumrows(env,lp);
        sprintf(cname[0], "flow(%d)", h);

        if (inst->turb_powers[h] < -0.5)
        {
            if (CPXnewrows(env, lp, 1, &zero, &less_equal, NULL, cname)) print_error("Wrong CPXnewrows() on flow constraints.");
            for (j = 0; j < inst->n_turbines; j++)
                if (CPXchgcoef(env, lp, lastrow, fpos(h, j, inst), 1.0)) print_error("Wrong CPXchgcoef() on flow constraints.");
            for (i = 0; i < inst->n_turbines; i++)
                if (CPXchgcoef(env, lp, lastrow, fpos(i, h, inst), -1.0)) print_error("Wrong CPXchgcoef() on flow constraints.");
        } else
        {
            if (CPXnewrows(env, lp, 1, &inst->turb_powers[h], &equal, NULL, cname)) print_error("Wrong CPXnewrows() on flow constraints.");
            for (j = 0; j < inst->n_turbines; j++)
                if (CPXchgcoef(env, lp, lastrow, fpos(h, j, inst), 1.0)) print_error("Wrong CPXchgcoef() on flow constraints.");
            for (i = 0; i < inst->n_turbines; i++)
                if (CPXchgcoef(env, lp, lastrow, fpos(i, h, inst), -1.0)) print_error("Wrong CPXchgcoef() on flow constraints.");
        }
    }

    // link constraints
    printf("%s Adding link constraints at row %d.\n", get_log(), CPXgetnumrows(env,lp));
    for (i = 0; i < inst->n_turbines; ++i)
    {
        for (j = 0; j < inst->n_turbines; j++)
        {
            lastrow = CPXgetnumrows(env,lp);
            //sprintf(cname[0], "link(%d,%d)", i+1, j+1);
            sprintf(cname[0], "link(%d,%d)", i, j);

            if (CPXnewrows(env, lp, 1, &zero, &equal, NULL, cname)) print_error("Wrong CPXnewrows() on link constraints.");
            for (k = 0; k < inst->n_cables; k++)
                if (CPXchgcoef(env, lp, lastrow, xpos(i, j, k, inst), 1.0)) print_error("Wrong CPXchgcoef() on link constraints.");

            if (CPXchgcoef(env, lp, lastrow, ypos(i, j, inst), -1.0)) print_error("Wrong CPXchgcoef() on link constraints.");
        }
    }

    // capacity link constraints
    printf("%s Adding capacity link constraints at row %d.\n", get_log(), CPXgetnumrows(env,lp));
    for (i = 0; i < inst->n_turbines; ++i)
    {
        for (j = 0; j < inst->n_turbines; j++)
        {
            lastrow = CPXgetnumrows(env,lp);
            //sprintf(cname[0], "capacity(%d,%d)", i+1, j+1);
            sprintf(cname[0], "capacity(%d,%d)", i, j);

            if (CPXnewrows(env, lp, 1, &zero, &greater_equal, NULL, cname)) print_error("Wrong CPXnewrows() on capacity link constraints.");
            for (k = 0; k < inst->n_cables; k++)
                if (CPXchgcoef(env, lp, lastrow, xpos(i, j, k, inst), inst->cable_powers[k])) print_error("Wrong CPXchgcoef() on capacity link constraints.");

            if (CPXchgcoef(env, lp, lastrow, fpos(i, j, inst), -1.0)) print_error("Wrong CPXchgcoef() on capacity link constraints.");
        }
    }

    // slack variables
    if (inst->slack_substation)
    {
        // substation in-degree constraint
        printf("%s Changing substation in-degree constraint at row %d.\n", get_log(), inst->indeg_constraint);
        if (CPXchgcoef(env, lp, inst->indeg_constraint, inst->sstart, -1.0)) print_error("Wrong CPXchgcoef() on slack variables.");
    }

    if (inst->slack_flow)
    {
        // flow constraints
        printf("%s Changing flow constraints at row %d.\n", get_log(), inst->flow_constraints_start);
        for (i = 0; i < inst->n_turbines; ++i)
        {
            if (CPXchgcoef(env, lp, inst->flow_constraints_start + i, inst->lstart + i, 1.0)) print_error("Wrong CPXchgcoef() on slack variables.");
        }
    }

    // get number of colons
    inst->ncols = CPXgetnumcols(env, lp);

    free(cname[0]);
    free(cname);
}

/*************************************************************************************************/
void build_model_cross_constraints(instance *inst, CPXENVptr env, CPXLPptr lp)
/*************************************************************************************************/
{
    // build the default model
    build_model(inst, env, lp);

    if (inst->use_cross_table)
    {
        inst->cross_table = generate_cross_table(inst);
    }

    CPXsetintparam(env,CPX_PARAM_REDUCE, CPX_PREREDUCE_NOPRIMALORDUAL);

    // add non-crossing contraints
    printf("%s Adding non-crossing contraints at row %d.\n", get_log(), CPXgetnumrows(env, lp));
    int ncuts = compute_no_cross_cuts(inst, env, lp, -1);
    printf("%s Added %d lazy (or new rows) constraints.\n", get_log(), ncuts);

    inst->ncols = CPXgetnumcols(env, lp);

    if (inst->use_cross_table)
    {
        free(inst->cross_table);
    }
}

/*************************************************************************************************/
void loop_solver(instance *inst, CPXENVptr env, CPXLPptr lp)
/*************************************************************************************************/
{
    int i, j, k, h, nzcnt, lastrow;
    int iteration = 0, newcuts = 0, ntotcuts = 0;
    double one = 1.0, time_diff;
    char less_equal = 'L';

    int *indeces = (int *) calloc(inst->ncols, sizeof(int));
    double *values = (double *) calloc(inst->ncols, sizeof(double));

    inst->model_type = LOOP_TASK_MODEL_TYPE;

    if (inst->use_cross_table)
    {
        inst->cross_table = generate_cross_table(inst);
    }
    
    while (1)
    {
        if (is_time_limit_expired(inst))
        {
            printf("%s LOOP TASK: time limit expired.\n", get_log());
            break;
        }

        time_diff = inst->time_limit - get_time_elapsed(inst);
        CPXsetdblparam(env, CPX_PARAM_TILIM, dmin(inst->iter_time, time_diff));		// try to solve it in 30 seconds
        CPXmipopt(env, lp);

        if (CPXgetbestobjval(env, lp, &inst->best_lb))
        {
            printf("%s LOOP TASK: no solution found.\n", get_log());
            continue;	// no solution found
        }

        CPXgetx(env, lp, inst->xstar, 0, inst->ncols - 1);
        newcuts = compute_no_cross_cuts(inst, env, lp, -1);
        ntotcuts += newcuts;

        if (inst->debug)
        {
            CPXwriteprob(env, lp, "model.lp", NULL);
        }

        clock_gettime(CLOCK_MONOTONIC, inst->time_end);

        //printf("\n");
        //printf("%s LOOP TASK\n\n", get_log());
        //printf("%s iteration =     %10.4d\n", get_log(), iteration);
        //printf("%s new cuts =      %10.4d\n", get_log(), newcuts);
        //printf("%s total cuts =    %10.4d\n", get_log(), ntotcuts);
        //printf("%s time elapsed =  %10.4lf\n", get_log(), ((double)inst->time_end->tv_sec + 1.0e-9*inst->time_end->tv_nsec) - ((double)inst->time_start->tv_sec + 1.0e-9*inst->time_start->tv_nsec));
        //printf("%s lower bound =   %10.4lf\n\n", get_log(), inst->best_lb);
        iteration++;
    }

    if (inst->use_cross_table)
    {
        free(inst->cross_table);
    }

    free(indeces);
    free(values);
}

/*************************************************************************************************/
void callback_solver(instance *inst, CPXENVptr env, CPXLPptr lp)
/*************************************************************************************************/
{
    // build the default model
    build_model(inst, env, lp);

    int i;

    CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF); 					                   // let MIP callbacks work on the original model

    if (inst->multi)
    {
        CPXsetlazyconstraintcallbackfunc(env, callback_cross_multi, (void *) inst);		       // install callback
    } else
    {
        CPXsetlazyconstraintcallbackfunc(env, callback_cross_single, (void *) inst);            // install callback        
    }
    CPXsetintparam(env, CPX_PARAM_THREADS, inst->num_threads); 					                // set core number

    if (inst->use_cross_table)
    {
        inst->cross_table = generate_cross_table(inst);
    }

    if (inst->multi)
    {
        inst->xstar_vector = (double **) malloc(inst->num_threads * sizeof(double *));
        for (i = 0; i < inst->num_threads; i++)
            inst->xstar_vector[i] = (double *) malloc(inst->ncols * sizeof(double));

        inst->indeces_vector = (int **) malloc(inst->num_threads * sizeof(int *));
        for (i = 0; i < inst->num_threads; i++)
            inst->indeces_vector[i] = (int *) malloc(inst->ncols * sizeof(int));

        inst->values_vector = (double **) malloc(inst->ncols * sizeof(double *));
        for (i = 0; i < inst->num_threads; i++)
            inst->values_vector[i] = (double *) malloc(inst->ncols * sizeof(double));
    }

    // Solve model
    CPXmipopt(env, lp);

    // Free all structures used in callbacks
    if (inst->use_cross_table)
    {
        free(inst->cross_table);
    }

    if (inst->multi)
    {
        for (i = 0; i < inst->num_threads; i++)
            free(inst->xstar_vector[i]);

        for (i = 0; i < inst->num_threads; i++)
            free(inst->indeces_vector[i]);

        for (i = 0; i < inst->num_threads; i++)
            free(inst->values_vector[i]);
    }
}

/*************************************************************************************************/
void hardfix_solver(instance *inst, CPXENVptr env, CPXLPptr lp)
/*************************************************************************************************/
{
    int count = 0;
    double prev_obj_val = DBL_MAX, time_diff;
    inst->obj_val = DBL_MAX;

    do
    {
        prev_obj_val = inst->obj_val;

        time_diff = inst->time_limit - get_time_elapsed(inst);
        CPXsetdblparam(env, CPX_PARAM_TILIM, dmin(inst->iter_time, time_diff));				// try to solve it in 30 seconds
        CPXmipopt(env, lp);

        CPXgetbestobjval(env, lp, &inst->best_lb);
        CPXgetobjval(env, lp, &inst->obj_val);
        
        CPXgetx(env, lp, inst->xstar, 0, inst->ncols-1);
        unfix_variables(inst, env, lp);
        hard_fix(inst, env, lp);
        //CPXwriteprob(env,lp, "model.lp", NULL);
        
        //printf("\n%s HARDFIX, iteration: %d\n", get_log(), count + 1);
        //printf("%s    Best integer=%lf\n", get_log(), inst->obj_val);
        //printf("%s    Best obj-val=%lf\n", get_log(), inst->best_lb);
        //printf("%s    Prev obj-val=%lf\n", get_log(), prev_obj_val);
        //printf("%s    Time elapsed=%lf\n", get_log(), get_time_elapsed(inst));
        
        count++;
    } while (prev_obj_val >= (inst->obj_val + EPSILON_SMALL) && !is_time_limit_expired(inst)); // TODO trovare condizione

    //CPXgetbestobjval(env, lp, &inst->best_lb);
    printf("\n%s Best obj-val=%lf\n", get_log(), inst->best_lb);
    //CPXgetobjval(env, lp, &inst->obj_val);
    printf("%s Best integer=%lf\n", get_log(), inst->obj_val);
}