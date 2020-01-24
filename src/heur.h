#ifndef HEUR_H_
#define HEUR_H_

#include "aux.h"
#include <math.h>

int dijkstra(instance *inst, int *preds);
int revisited_dijkstra(instance *inst, int *preds);
int preds_to_solution(instance *inst, int *preds, double *solution);
int preds_to_solution_params(instance *inst, int *preds, double *solution, int i, int j);
double VNS(instance *inst, int *preds, double best_score);
double simulated_annealing(instance *inst, int *preds, double best_score);
double one_opt_move(instance *inst, int *preds, double best_score);
double two_opt_move(instance *inst, int *preds, double best_score);
double three_opt_move(instance *inst, int *preds, double best_score);
double evaluate_preds(instance *inst, int *preds, double best_score, int log);
double evaluate_preds_param(instance *inst, int *preds, int a, int b, int c, double best_score, int log);
double evaluate_solution(instance *inst, double *solution, double best_score, int log);
int heur_VNS_launcher(instance *inst);
int heur_simulated_annealing_launcher(instance *inst);

#endif // HEUR_H_
