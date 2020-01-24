
#include "wf.h"

void main_free_instance(instance *inst);
void main_parse_command_line(int argc, char** argv, instance *inst);
void main_read_input(instance *inst);

instance inst;

int main(int argc, char **argv) {

    // No enought arguments provided
    if (argc < 2) {
        printf("Usage: %s -help for help\n", argv[0]);
        exit(1);
    }

    main_parse_command_line(argc, argv, &inst);
    main_read_input(&inst);

    struct timespec tstart, tend;
    inst.time_start = &tstart;
    inst.time_end = &tend;

    clock_gettime(CLOCK_MONOTONIC, inst.time_start);
    if (WFopt(&inst)) print_error("Error within WFopt().");
    clock_gettime(CLOCK_MONOTONIC, inst.time_end);

    if (VERBOSE >= 1) printf("%s WFopt ended. Time elapsed = %lf seconds.\n", get_log(), ((double)inst.time_end->tv_sec + 1.0e-9*inst.time_end->tv_nsec) - ((double)inst.time_start->tv_sec + 1.0e-9*inst.time_start->tv_nsec));

    main_free_instance(&inst);
    return 0;
}

void main_free_instance(instance *inst) {
    free(inst->x_turb_coords);
    free(inst->y_turb_coords);
    free(inst->turb_powers);
    free(inst->cable_powers);
    free(inst->cable_costs);
}

void main_parse_command_line(int argc, char** argv, instance *inst) {
    int help, i;
    help = 0;

    // default
    strcpy(inst->turbine_input_file, "NULL");
    strcpy(inst->cable_input_file, "NULL");
    strcpy(inst->params_file, "NULL");
    inst->subs_cables = 999999999;
    inst->turb_tot_power = 0.0;

    inst->slack_substation = 0;
    inst->slack_flow = 0;
    inst->print_model = 0;
    inst->rins = 0;
    inst->polishafter=1e+75; 					// default value for polishafter
    strcpy(inst->solver, "NULL");
    inst->model_type = 0;
    inst->polishafter = 0.0;
    inst->randomseed = 0;
    inst->use_cross_table = 0;
    inst->multi = 0;
    inst->num_threads = 4;
    inst->time_limit = CPX_INFBOUND;
    inst->iter_time = 30.0;
    inst->hard_fix_prob = 0.5;

    inst->time_limit_expired = 0;
    
    if (argc < 1) help = 1;

    for (i = 1; i < argc; i++)
    {
        if (strcmp(argv[i],"-turbine_file") == 0) { strcpy(inst->turbine_input_file, argv[++i]); continue; }
        if (strcmp(argv[i],"-cable_file") == 0) { strcpy(inst->cable_input_file, argv[++i]); continue; }
        if (strcmp(argv[i],"-subs_cables") == 0) { inst->subs_cables = atoi(argv[++i]); continue; }
        if (strcmp(argv[i],"-print_model") == 0) { inst->print_model = 1; continue; }
        if (strcmp(argv[i],"-rins") == 0) { inst->rins = atoi(argv[++i]); continue; }						// rins
        if (strcmp(argv[i],"-polishafter") == 0) { inst->polishafter = atof(argv[++i]); continue; }			// polishafter
        if (strcmp(argv[i],"-time_limit") == 0) { inst->time_limit = atof(argv[++i]); continue; }           // total time limit
        if (strcmp(argv[i],"-iter_time") == 0) { inst->iter_time = atof(argv[++i]); continue; }             // iteration time limit
        if (strcmp(argv[i],"-hard_fix_prob") == 0) { inst->hard_fix_prob = atof(argv[++i]); continue; }     // hard fix probability
        if (strcmp(argv[i],"-solver") == 0) { strcpy(inst->solver, argv[++i]); continue; }
        if (strcmp(argv[i],"-slack_substation") == 0) { inst->slack_substation = 1; continue; } 			// add slack variable substation constraint
        if (strcmp(argv[i],"-slack_flow") == 0) { inst->slack_flow = 1; continue; }                         // add slack variables flow constraints
        if (strcmp(argv[i],"-seed") == 0) { inst->randomseed = abs(atoi(argv[++i])); continue; } 			// random seed
        if (strcmp(argv[i],"-use_cross_table") == 0) { inst->use_cross_table = 1; continue; }               // use cross talbe
        if (strcmp(argv[i],"-multi") == 0) { inst->multi = 1; continue; }                                   // use multi-thread code
        if (strcmp(argv[i],"-num_threads") == 0) { inst->num_threads = atoi(argv[++i]); continue; } 		// n. threads
        if (strcmp(argv[i],"-debug") == 0) { inst->debug = 1; continue; }
        if (strcmp(argv[i],"-help") == 0) { help = 1; continue; }											// help
    }

    if (help || (VERBOSE >= 10))				// print current parameters
    {
        printf("\n*--- Available parameters (vers. 23-mar-2018) ---------------*\n\n");
        printf("        turbine_file:   %s\n", inst->turbine_input_file);
        printf("        cable_file:     %s\n", inst->cable_input_file);
        printf("        subs_cables:    %d\n", inst->subs_cables);
        printf("        print_model:    %d\n", inst->print_model);
        printf("        rins:           %d\n", inst->rins);
        printf("        polishafter:    %lf\n", inst->polishafter);
        printf("        time_limit:     %lf\n", inst->time_limit);
        printf("        iter_time:      %lf\n", inst->iter_time);
        printf("        hard_fix_prob:  %lf\n", inst->hard_fix_prob);
        printf("        solver:         %s\n", inst->solver);
        printf("        slack_subs:     %d\n", inst->slack_substation);
        printf("        slack_flow:     %d\n", inst->slack_flow);
        printf("        seed:           %d\n", inst->randomseed);
        printf("        cross_table:    %d\n", inst->use_cross_table);
        printf("        multi:          %d\n", inst->multi);
        printf("        threads:        %d\n", inst->num_threads);
        printf("        enter -help for help\n\n");
        printf("*--------------------------------------------------------------*\n\n");
    }

    if (help) exit(1);
}

void main_read_input(instance *inst) {
    char line[100];
    int index, line_counter;
    int x, y, p;
    double c, l;

    FILE *fin;

    // TURBINES
    fin = fopen(inst->turbine_input_file, "r");
    if (fin == NULL) print_error("Turbine input file not found.");

    line_counter = 0;

    // get number of turbine entries
    while (fgets(line, sizeof(line), fin) != NULL)
    {
        if (feof(fin)) break;
        line_counter++;
    }
    inst->n_turbines = line_counter;
    fclose(fin);

    inst->x_turb_coords = (double *) calloc(inst->n_turbines, sizeof(double));
    inst->y_turb_coords = (double *) calloc(inst->n_turbines, sizeof(double));
    inst->turb_powers = (double *) calloc(inst->n_turbines, sizeof(double));

    // read the file and insert data
    fin = fopen(inst->turbine_input_file, "r");
    for (index = 0; index < inst->n_turbines; index++)
    {
        if (!fscanf(fin, "%d %d %d", &x, &y, &p)) print_error("Wrong turbine file format.");

        inst->x_turb_coords[index] = (double) x;
        inst->y_turb_coords[index] = (double) y;
        inst->turb_powers[index] = (double) p;

        // looking for the substation, where power is -1
        if (p < -0.5)
        {
            inst->subs_index = index;
            inst->x_subs_coord = (double) x;
            inst->y_subs_coord = (double) y;
        } else
        {
            inst->turb_tot_power += p;
        }

        if (VERBOSE >= 1000) printf("x=%lf, y=%lf, power=%lf\n", inst->x_turb_coords[index], inst->y_turb_coords[index], inst->turb_powers[index]);
    }
    fclose(fin);


    // CABLES
    fin = fopen(inst->cable_input_file, "r");
    if (fin == NULL) print_error("Cable input file not found.");

    line_counter = 0;

    // get number of cable entries
    while (fgets(line, sizeof(line), fin) != NULL)
    {
        if (feof(fin)) break;
        line_counter++;
    }
    inst->n_cables = line_counter;
    fclose(fin);

    inst->cable_powers = (double *) calloc(inst->n_cables, sizeof(double));
    inst->cable_costs = (double *) calloc(inst->n_cables, sizeof(double));

    // read the file and insert data
    fin = fopen(inst->cable_input_file, "r");
    for (index = 0; index < inst->n_cables; index++)
    {
        if (!fscanf(fin, "%d %lf %lf", &p, &c, &l)) print_error("Wrong cable file format.");

        inst->cable_powers[index] = p;
        inst->cable_costs[index] = c;

        if ( VERBOSE >= 1000 ) printf("power=%lf, cost=%lf\n", inst->cable_powers[index], inst->cable_costs[index]);
    }
    fclose(fin);
}