
#include "wf_heur.hpp"

int main_parse_command_line(int argc, char **argv, instance *inst);

instance inst;

int main(int argc, char *argv[])
{
    // No enought arguments provided
    if (argc < 2)
    {
        printf("Usage: %s --help to display options\n", argv[0]);
        exit(1);
    }

    if (main_parse_command_line(argc, argv, &inst))
    {
        exit(1);
    }
    if (comm_read_input(&inst))
    {
        exit(1);
    }

    struct timespec tstart, tend;
    inst.time_start = &tstart;
    inst.time_end = &tend;

    clock_gettime(CLOCK_MONOTONIC, inst.time_start);
    if (inst.solver.compare("vns") == 0)
    {
        wf_vns_launcher(&inst);
    }
    else if (inst.solver.compare("siman") == 0)
    {
        wf_siman_launcher(&inst);
    }
    else
    {
        comm_print_error("Solver unknown.");
    }
    clock_gettime(CLOCK_MONOTONIC, inst.time_end);

    if (VERBOSE >= 1)
        printf("WFopt ended. Time elapsed = %lf seconds.\n", ((double)inst.time_end->tv_sec + 1.0e-9 * inst.time_end->tv_nsec) - ((double)inst.time_start->tv_sec + 1.0e-9 * inst.time_start->tv_nsec));

    comm_free_instance(&inst);
    return 0;
}

int main_parse_command_line(int argc, char *argv[], instance *inst)
{
    try
    {
        po::options_description desc("Allowed options");
        desc.add_options()
        ("help", "produce help message")
        ("turbines", po::value<string>(&inst->turbine_input_file), "set turbine file path")
        ("cables", po::value<string>(&inst->cable_input_file), "set cable file path")
        ("subs_cables", po::value<int>(&inst->subs_cables)->default_value(999999), "set substation input cables number")
        ("solver", po::value<string>(&inst->solver)->default_value("none"), "set solver: 'vns', 'siman'")
        //("print_model", po::value<bool>(&inst->print_model)->default_value(false), "print LP model after built")
        ("time_limit", po::value<double>(&inst->time_limit)->default_value(3600), "set solver time limit (s)")
        //("iter_time", po::value<double>(&inst->iter_time)->default_value(30.0), "set iteration time (used for ...)")
        //("slack_subs", po::value<bool>(&inst->slack_substation)->default_value(false), "set slack for substation")
        //("slack_flow", po::value<bool>(&inst->slack_flow)->default_value(false), "set slack variables for flows")
        ("seed", po::value<int>(&inst->randomseed)->default_value(0), "set random seed")
        ("cross_table", po::value<bool>(&inst->use_cross_table)->default_value(false), "set cross table optimization")
        ("multi", po::value<bool>(&inst->multi)->default_value(false), "set multi-thread code")
        ("threads", po::value<int>(&inst->num_threads)->default_value(4), "set number of threads");

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if (vm.count("help"))
        {
            cout << desc << "\n";
            return 1;
        }

        if (vm.count("turbines"))
        {
            cout << "Turbine file path set to " << vm["turbines"].as<string>() << ".\n";
        }
        else
        {
            cout << "Turbine file path not set. Exiting.\n";
            return 1;
        }

        if (vm.count("cables"))
        {
            cout << "Cable file path set to " << vm["cables"].as<string>() << ".\n";
        }
        else
        {
            cout << "Cable file path not set. Exiting.\n";
            return 1;
        }

        if (vm.count("subs_cables"))
        {
            cout << "Subs cables set to " << vm["subs_cables"].as<int>() << ".\n";
        }
        if (vm.count("solver"))
        {
            cout << "Solver set to " << vm["solver"].as<string>() << ".\n";
        }
        if (vm.count("time_limit"))
        {
            cout << "Time limit set to " << vm["time_limit"].as<double>() << ".\n";
        }
        if (vm.count("iter_time"))
        {
            cout << "Iteration time set to " << vm["iter_time"].as<double>() << ".\n";
        }
        if (vm.count("print_model"))
        {
            cout << "Print model set to " << vm["print_model"].as<bool>() << ".\n";
        }
        if (vm.count("slack_subs"))
        {
            cout << "Slack subs set to " << vm["slack_subs"].as<bool>() << ".\n";
        }
        if (vm.count("slack_flow"))
        {
            cout << "Slack flow set to " << vm["slack_flow"].as<bool>() << ".\n";
        }
        if (vm.count("cross_table"))
        {
            cout << "Cross table set to " << vm["cross_table"].as<bool>() << ".\n";
        }
        if (vm.count("multi"))
        {
            cout << "Multi set to " << vm["multi"].as<bool>() << ".\n";
        }
        if (vm.count("seed"))
        {
            cout << "Random seed set to " << vm["seed"].as<int>() << ".\n";
        }
        if (vm.count("threads"))
        {
            cout << "Number of threads set to " << vm["threads"].as<int>() << ".\n";
        }
    }
    catch (exception &e)
    {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch (...)
    {
        cerr << "Exception of unknown type!\n";
        return 1;
    }

    return 0;
}
