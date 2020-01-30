# Wind Farm Cables Routing Optimization

The software consists of a solver for finding a minimum cost interconnection between all the turbines of a given wind farm, respecting the specific design constraints.

This project implements several approaches that work with the **Mixed-Integer Linear Programming** (MILP) solver CPLEX.
Alternatively, we implement two heuristic solvers that can run independently from the CPLEX library.

The original article we utilised to set up this project can be found [here](https://orbit.dtu.dk/en/publications/optimizing-wind-farm-cable-routing-considering-power-losses-2).
The results obtained are described in this [report](https://github.com/caerbannogwhite/WindFarmCablesRoutingOptimization/blob/master/report/report.pdf).

The software requires [CPLEX](https://www.ibm.com/analytics/cplex-optimizer) and the [Boost](https://www.boost.org/) library. Once you have installed them in your computer, you should change the `CPLEX_HOME` and the `BOOST_HOME` variables in `Makefile`.
If `Makefile` is not able to locate `CPLEX_HOME` correctly, then it will produce only one executable file (`wfsolver_heur`). Otherwise, it will also produce the `wfsolver_cpx` executable.

Here are some images of the connections computed by our solver.

[Solution for testbed 03 obtained with the *hardfix solver*.](https://github.com/caerbannogwhite/WindFarmCablesRoutingOptimization/blob/master/report/img/hardfix_03_p03_res.png)


[Solution for testbed 27 obtained with the *hardfix solver*.](https://github.com/caerbannogwhite/WindFarmCablesRoutingOptimization/blob/master/report/img/hardfix_27_p03_res.png)

The original testbeds are reported in the following table.

|         TURBINES |            CABLES | SUBSTATION
|------------------|-------------------|-----------
| data/data_01.cbl | data/data_01.turb | 10
| data/data_02.cbl | data/data_02.turb | 10
| data/data_03.cbl | data/data_03.turb | 10
| data/data_04.cbl | data/data_04.turb | 10
| data/data_05.cbl | data/data_05.turb | 10
| data/data_06.cbl | data/data_06.turb | 10
| data/data_07.cbl | data/data_07.turb | 1000
| data/data_08.cbl | data/data_08.turb | 1000
| data/data_09.cbl | data/data_09.turb | 1000
| data/data_10.cbl | data/data_10.turb | 1000
| data/data_12.cbl | data/data_12.turb | 1000
| data/data_13.cbl | data/data_13.turb | 1000
| data/data_14.cbl | data/data_14.turb | 1000
| data/data_15.cbl | data/data_15.turb | 1000
| data/data_16.cbl | data/data_16.turb | 4
| data/data_17.cbl | data/data_17.turb | 4
| data/data_18.cbl | data/data_18.turb | 4
| data/data_19.cbl | data/data_19.turb | 4
| data/data_20.cbl | data/data_20.turb | 10
| data/data_21.cbl | data/data_21.turb | 10
| data/data_26.cbl | data/data_26.turb | 10
| data/data_27.cbl | data/data_27.turb | 10
| data/data_28.cbl | data/data_28.turb | 10
| data/data_29.cbl | data/data_29.turb | 10
