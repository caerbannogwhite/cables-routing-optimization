
#ifndef WF_HEUR_H_
#define WF_HEUR_H_

#include "common.hpp"

int dijkstra(instance *inst, int *preds);
int revisited_dijkstra(instance *inst, int *preds);
int preds_to_solution(instance *inst, int *preds, double *solution);
int preds_to_solution_params(instance *inst, int *preds, double *solution, int i, int j);
double wf_vns(instance *inst, int *preds, double best_score);
double wf_siman(instance *inst, int *preds, double best_score);
double one_opt_move(instance *inst, int *preds, double best_score);
double two_opt_move(instance *inst, int *preds, double best_score);
double three_opt_move(instance *inst, int *preds, double best_score);
double evaluate_preds(instance *inst, int *preds, double best_score, int log);
double evaluate_preds_param(instance *inst, int *preds, int a, int b, int c, double best_score, int log);
double evaluate_solution(instance *inst, double *solution, double best_score, int log);
int wf_vns_launcher(instance *inst);
int wf_siman_launcher(instance *inst);

#endif // WF_HEUR_H_
