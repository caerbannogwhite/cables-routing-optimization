
#include "wf.h"

int main_free_instance(instance *inst);
int main_parse_command_line(int argc, char** argv, instance *inst);
int main_read_input(instance *inst);

instance inst;

int main(int argc, char *argv[])
{
    // No enought arguments provided
    if (argc < 2) {
        printf("Usage: %s --help to display options\n", argv[0]);
        exit(1);
    }

    if (main_parse_command_line(argc, argv, &inst)) { exit(1); }
    if (main_read_input(&inst)) { exit(1); }

    struct timespec tstart, tend;
    inst.time_start = &tstart;
    inst.time_end = &tend;

    clock_gettime(CLOCK_MONOTONIC, inst.time_start);
    if (WFopt(&inst)) print_error("Error within WFopt().");
    clock_gettime(CLOCK_MONOTONIC, inst.time_end);

    if (VERBOSE >= 1) printf("WFopt ended. Time elapsed = %lf seconds.\n", ((double)inst.time_end->tv_sec + 1.0e-9*inst.time_end->tv_nsec) - ((double)inst.time_start->tv_sec + 1.0e-9*inst.time_start->tv_nsec));

    main_free_instance(&inst);
    return 0;
}

int main_free_instance(instance *inst)
{
    free(inst->x_turb_coords);
    free(inst->y_turb_coords);
    free(inst->turb_powers);
    free(inst->cable_powers);
    free(inst->cable_costs);

    return 0;
}

int main_parse_command_line(int argc, char* argv[], instance *inst)
{
    try
    {
        po::options_description desc("Allowed options");
        desc.add_options()
        ("help", "produce help message")
        ("turbines", po::value<string>(&inst->turbine_input_file), "set turbine file path")
        ("cables", po::value<string>(&inst->cable_input_file), "set cable file path")
        ("subs_cables", po::value<int>(&inst->subs_cables)->default_value(999999), "set substation input cables number")
        ("solver", po::value<string>(&inst->solver)->default_value("none"), "set solver")
        ("print_model", po::value<bool>(&inst->print_model)->default_value(false), "print LP model after built")
        ("time_limit", po::value<double>(&inst->time_limit)->default_value(CPX_INFBOUND), "set solver time limit")
        ("iter_time", po::value<double>(&inst->iter_time)->default_value(30.0), "set iteration time (used for ...)")
        ("hard_fix_prob", po::value<double>(&inst->hard_fix_prob)->default_value(0.5), "set hard fix probability")
        ("polishafter", po::value<double>(&inst->polishafter)->default_value(1e75), "set cplex polishafter value")
        ("rins", po::value<bool>(&inst->rins)->default_value(false), "if 1, set cplex rins value")
        ("slack_subs", po::value<bool>(&inst->slack_substation)->default_value(false), "set slack for substation")
        ("slack_flow", po::value<bool>(&inst->slack_flow)->default_value(false), "set slack variables for flows")
        ("seed", po::value<int>(&inst->randomseed)->default_value(0), "set random seed")
        ("cross_table", po::value<bool>(&inst->use_cross_table)->default_value(false), "set cross table optimization")
        ("multi", po::value<bool>(&inst->multi)->default_value(false), "set multi-thread code")
        ("threads", po::value<int>(&inst->num_threads)->default_value(4), "set number of threads");

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if (vm.count("help")) { cout << desc << "\n"; return 1; }

        if (vm.count("turbines")) { cout << "Turbine file path set to " << vm["turbines"].as<string>() << ".\n"; }
        else { cout << "Turbine file path not set. Exiting.\n"; return 1; }

        if (vm.count("cables")) { cout << "Cable file path set to " << vm["cables"].as<string>() << ".\n"; }
        else { cout << "Cable file path not set. Exiting.\n"; return 1; }

        if (vm.count("subs_cables")) { cout << "Subs cables set to " << vm["subs_cables"].as<int>() << ".\n"; }
        if (vm.count("solver")) { cout << "Solver set to " << vm["solver"].as<string>() << ".\n"; }
        if (vm.count("time_limit")) { cout << "Time limit set to " << vm["time_limit"].as<double>() << ".\n"; }
        if (vm.count("iter_time")) { cout << "Iteration time set to " << vm["iter_time"].as<double>() << ".\n"; }
        if (vm.count("hard_fix_prob")) { cout << "Hard fix probability set to " << vm["hard_fix_prob"].as<double>() << ".\n"; }
        if (vm.count("polishafter")) { cout << "CPLEX polishafter set to " << vm["polishafter"].as<double>() << ".\n"; }
        
        if (vm.count("print_model")) { cout << "Print model set to " << vm["print_model"].as<bool>() << ".\n"; }
        if (vm.count("rins")) { cout << "CPLEX rins set to " << vm["rins"].as<bool>() << ".\n"; }
        if (vm.count("slack_subs")) { cout << "Slack subs set to " << vm["slack_subs"].as<bool>() << ".\n"; }
        if (vm.count("slack_flow")) { cout << "Slack flow set to " << vm["slack_flow"].as<bool>() << ".\n"; }
        if (vm.count("cross_table")) { cout << "Cross table set to " << vm["cross_table"].as<bool>() << ".\n"; }
        if (vm.count("multi")) { cout << "Multi set to " << vm["multi"].as<bool>() << ".\n"; }

        if (vm.count("seed")) { cout << "Random seed set to " << vm["seed"].as<int>() << ".\n"; }
        if (vm.count("threads")) { cout << "Number of threads set to " << vm["threads"].as<int>() << ".\n"; }
    }
    catch (exception &e)
    {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch (...)
    {
        cerr << "Exception of unknown type!\n";
    }

    return 0;
}

int main_read_input(instance *inst)
{
    char line[100];
    int index, line_counter;
    int x, y, p;
    double c, l;

    FILE *fin;

    // TURBINES
    fin = fopen(&inst->turbine_input_file[0], "r");
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
    fin = fopen(&inst->turbine_input_file[0], "r");
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
    fin = fopen(&inst->cable_input_file[0], "r");
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
    fin = fopen(&inst->cable_input_file[0], "r");
    for (index = 0; index < inst->n_cables; index++)
    {
        if (!fscanf(fin, "%d %lf %lf", &p, &c, &l)) print_error("Wrong cable file format.");

        inst->cable_powers[index] = p;
        inst->cable_costs[index] = c;

        if ( VERBOSE >= 1000 ) printf("power=%lf, cost=%lf\n", inst->cable_powers[index], inst->cable_costs[index]);
    }
    fclose(fin);

    return 0;
}