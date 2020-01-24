# Wind Farm Cables Routing Optimization

We developed this project for the Operations Research 2 course.
The software consists of a solver for finding a minimum cost interconnection between all the turbines of a given wind farm, respecting the specific design constraints.

The software implements several approaches that work with the Mixed Integer Linear Programming (MILP) solver CPLEX.
Alternatively, we implement two heuristic solvers that can run independently from the CPLEX library.

The original article we utilised to set up this project can be found [here](https://orbit.dtu.dk/en/publications/optimizing-wind-farm-cable-routing-considering-power-losses-2).

The software requires **CPLEX**. Once you have installed it in your computer, you can change the `CPLEX_HOME` variable in `Makefile`.

Here are some images of the connections computed by our solver.

[hardfix_03](https://github.com/caerbannogwhite/WindFarmCablesRoutingOptimization/blob/master/report/img/hardfix_03_p03_res.png)
[hardfix_27](https://github.com/caerbannogwhite/WindFarmCablesRoutingOptimization/blob/master/report/img/hardfix_27_p03_res.png)