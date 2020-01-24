
#include "heur.h"

int dijkstra(instance *inst, int *preds) {
    int i, j, k, h, sign;
    double c, dist;

    int *flags = (int *) calloc(sizeof(int), inst->n_turbines);
    double *L = (double *) calloc(sizeof(double), inst->n_turbines);
    double *costs = (double *) calloc(sizeof(double), inst->n_turbines * inst->n_turbines);
    
    double perturbation = 0.99;
    int start = inst->subs_index;

    // initialize perturbated costs
    h = 0;
    for (i = 0; i < inst->n_turbines; ++i)
    {
        for (j = i + 1; j < inst->n_turbines; ++j)
        {
            // perturbation
            sign = (rand() + 0.0) / RAND_MAX > 0.5 ? 1 : -1;
            dist = hypot(inst->x_turb_coords[i] - inst->x_turb_coords[j], inst->y_turb_coords[i] - inst->y_turb_coords[j]) * (1 + perturbation * sign);
            
            // completely random
            //dist = rand() % 80 + 20;

            costs[i*inst->n_turbines + j] = dist;
            costs[j*inst->n_turbines + i] = dist;
        }
    }

    // dijkstra initialization
    for (j = 0; j < inst->n_turbines; ++j)
    {
        if(j != start)
        {
            flags[j] = 0;
            L[j] = costs[j*inst->n_turbines + start];
            preds[j] = start;
        }
        else 
        {
            flags[start] = 1;
            preds[start] = start;
            L[start] = 0.0;
        }
    }
    
    // dijkstra edge selection
    for (k = 0; k < inst->n_turbines; ++k)
    {
        double min = DBL_MAX;
        for (j = 0; j < inst->n_turbines; ++j) // select n-1 edges
        {
            // TODO: randomize!!!
            if((flags[j] == 0) && (L[j] < min))
            {
                min = L[j];
                h = j;
            }
        }

        flags[h] = 1; // extends S
        
        for (j = 0; j < inst->n_turbines; ++j) // update L[j] and pred[j] for every j not in S
        {
            c = costs[h*inst->n_turbines + j];
            if ((flags[j] == 0) && (c < L[j])) 
            {
                L[j] = c;
                preds[j] = h;
            }
        }
    }

    // TODO: aggiungere archi entranti nella substation

    free(flags);
    free(L);
    free(costs);

    return 0;
}


// returns 1 there is a cycle
int preds_to_solution(instance *inst, int *preds, double *solution) {
    int i, j, counter, edges;
    
    // clear and initialize structures
    // streams[i] is the input stream in the i-th turbine
    for (i = 0; i < inst->ncols; ++i) solution[i] = 0.0;

    for (i = 0; i < inst->n_turbines; ++i)
    {
        inst->streams[i] = 0.0;
        inst->preds_copy[i] = preds[i];
        inst->flags[i] = 1;
    }

    edges = 0;
    while (1)
    {
        for (i = 0; i < inst->n_turbines; ++i)
        {
            if (inst->preds_copy[i] == -1) inst->flags[i] = 0;
            inst->flags[inst->preds_copy[i]] = 0;
        }

        counter = 0;
        for (i = 0; i < inst->n_turbines; ++i)
        {
            if (inst->flags[i])
            {
                inst->streams[inst->preds_copy[i]] += inst->streams[i] + 1.0;
                inst->preds_copy[i] = -1;
                ++counter;
            }
        }

        // ends correctly
        if (counter == 0) break;

        // there is a cycle
        if (edges > inst->n_turbines)
        {
            free(inst->preds_copy);
            free(inst->flags);
            free(inst->streams);
            return 1;
        }
        ++edges;

        for (i = 0; i < inst->n_turbines; ++i) inst->flags[i] = 1;
    }

    for (i = 0; i < inst->n_turbines; ++i)
    {
        if (i == preds[i]) continue;

        solution[fpos(i, preds[i], inst)] = inst->streams[i] + 1;
        solution[ypos(i, preds[i], inst)] = 1.0;

        for (j = 0; j < inst->n_cables; ++j)
        {
            if (inst->cable_powers[j] >= inst->streams[i] + 1) break;
        }

        if (inst->cable_powers[inst->n_cables - 1] < inst->streams[i] + 1)
        {
            solution[xpos(i, preds[i], inst->n_cables - 1, inst)] = 1.0;
        }
        else
        {
            solution[xpos(i, preds[i], j, inst)] = 1.0; 
        }
    }

    return 0;
}


double VNS(instance *inst, int *preds, double best_score) {
    int i, oom_improvement = 0;
    double oom_best_score, eval_preds;

    //best_score = evaluate_preds(inst, preds, DBL_MAX, 0);
    //printf("start: best_score=%lf\n", best_score);

    while (1)
    {
        oom_best_score = one_opt_move(inst, preds, best_score);

        // Local optimum not reached yet
        if (best_score > oom_best_score)
        {
            //printf("%s     Optimum not reaced: best_score = %.4e, oom_score = %.4e\n", get_log(), best_score, oom_best_score);
            best_score = oom_best_score;
            oom_improvement = 1;

            //printf("%s time_elapsed=%lf Solution found with score=%.4e\n", get_log(), get_time_elapsed(inst), best_score);
        }

        // Local optimum reached: launching 3-opt move
        else
        {
            //printf("%s   Local optimum found with score = %.4e. Launching TOM.\n", get_log(), best_score);
            if (oom_improvement) {
                //score = three_opt_move(inst, preds, best_score);
                oom_improvement = 0;
            }
            else {
                break;
            }
            
            // DEBUG
            eval_preds = evaluate_preds(inst, preds, DBL_MAX, 0);
            //double *sol = (double *) calloc(inst->ncols, sizeof(double));
            //preds_to_solution(inst, preds, sol);
            //double eval_sol = evaluate_solution(inst, inst->solution, DBL_MAX, 0);
            //printf("tom_score=%lf, eval_preds=%lf, eval_sol=%lf\n",score, eval_preds, eval_sol);
            //free(sol);
            if (best_score > eval_preds) {
                for (i = 0; i < inst->ncols; ++i) inst->xstar[i] = inst->solution[i];
                best_score = eval_preds;

                printf("%s time_elapsed=%lf Solution found with score=%.4e\n", get_log(), get_time_elapsed(inst), best_score);
            }
        }

        if (inst->time_limit_expired) break;
    }
    return best_score;
}


double one_opt_move(instance *inst, int *preds, double best_score) {
    int i, j, tmp1;
    double score;
    
    for (i = 0; i < inst->n_turbines; ++i)
    {
        for (j = 0; j < inst->n_turbines; ++j)
        {
            tmp1 = preds[i];
            preds[i] = j;

            score = evaluate_preds(inst, preds, best_score, 0);
                                    
            if (best_score > score)
            {
                return score;
            } else {
                preds[i] = tmp1;
            }
            if (is_time_limit_expired(inst)) break;

            if (inst->time_limit_expired) break;
        }

        if (inst->time_limit_expired) break;
    }

    return best_score;
}


double two_opt_move(instance *inst, int *preds, double best_score) {
	int i, j, k, tmp1, tmp2;
    double score;
    
    for (i = 0; i < inst->n_turbines; ++i)
    {
        for (j = 0; j < inst->n_turbines; ++j)
        {
            for (k = 0; k < inst->n_turbines; ++k)
            {
                tmp1 = preds[i];
                tmp2 = preds[j];

                preds[i] = j;
                preds[j] = k;

                score = evaluate_preds(inst, preds, best_score, 0);
                                        
                if ((best_score - score) > EPSILON_MED)
                {
                    return score;
                } else {
                    preds[i] = tmp1;
                    preds[j] = tmp2;
                }
                if (is_time_limit_expired(inst)) break;
            }

            if (inst->time_limit_expired) break;
        }

        if (inst->time_limit_expired) break;
    }

    return best_score;
}


double three_opt_move(instance *inst, int *preds, double best_score) {
    int i, j, k, h, tmp1, tmp2;
    double score;

    for (i = 0; i < inst->n_turbines; ++i)
    {
        for (j = 0; j < inst->n_turbines; ++j)
        {
            for (k = 0; k < inst->n_turbines; ++k)
            {
                for (h = 0; h < inst->n_turbines; ++h)
                {
                    tmp1 = preds[i];
                    tmp2 = preds[k];

                    preds[i] = j;
                    preds[k] = h;

                    score = evaluate_preds(inst, preds, best_score, 0);
                                        
                    if (best_score > score)
                    {
                        return score;
                    } else {
                        preds[i] = tmp1;
                        preds[k] = tmp2;
                    }
                }

                if (is_time_limit_expired(inst)) break;
            }

            if (inst->time_limit_expired) break;
        }

        if (inst->time_limit_expired) break;
    }

    return best_score;
}


/*  evaluate_preds: valuta una soluzione nella forma della
    struttura a predecessori, fermandosi non appena la
    valutazione oltrepassa la soglia individuata dal parametro
    best_score.

    Strutture utilizzate:

    inst->flags
    inst->preds_copy
    inst->solution
    inst->streams
*/
double evaluate_preds(instance *inst, int *preds, double best_score, int log) {
    int cnt, edges;
    int i, j, k, h, in_deg_cnt = 0, cross_cnt = 0, stream_exceed_cnt = 0;
    double score = 0.0;
    
    // clear and initialize structures
    // streams[i] is the input stream in the i-th turbine
    for (i = 0; i < inst->ncols; ++i) inst->solution[i] = 0.0;

    for (i = 0; i < inst->n_turbines; ++i)
    {
        inst->streams[i] = 0.0;
        inst->preds_copy[i] = preds[i];
        inst->flags[i] = 1;
    }

    edges = 0;
    while (1)
    {
        for (i = 0; i < inst->n_turbines; ++i)
        {
            if (inst->preds_copy[i] == -1) inst->flags[i] = 0;
            inst->flags[inst->preds_copy[i]] = 0;
        }

        cnt = 0;
        for (i = 0; i < inst->n_turbines; ++i)
        {
            if (inst->flags[i])
            {
                inst->streams[inst->preds_copy[i]] += inst->streams[i] + 1.0;
                inst->preds_copy[i] = -1;
                ++cnt;
            }
        }

        // ends correctly
        if (cnt == 0) break;

        // there is a cycle
        if (edges > inst->n_turbines)
        {
            return DBL_MAX;
        }
        ++edges;

        for (i = 0; i < inst->n_turbines; ++i) inst->flags[i] = 1;
    }

    if ((inst->streams[0] - inst->turb_tot_power) > EPSILON_SMALL || (inst->streams[0] - inst->turb_tot_power) < -EPSILON_SMALL)
    {
        // not acceptable
        return DBL_MAX;
    }

    for (i = 0; i < inst->n_turbines; ++i)
    {
        if (i == preds[i]) continue;

        inst->solution[fpos(i, preds[i], inst)] = inst->streams[i] + 1;
        inst->solution[ypos(i, preds[i], inst)] = 1.0;

        for (j = 0; j < inst->n_cables; ++j)
        {
            if (inst->cable_powers[j] >= inst->streams[i] + 1) break;
        }

        if (inst->cable_powers[inst->n_cables - 1] < inst->streams[i] + 1)
        {
            inst->solution[xpos(i, preds[i], inst->n_cables - 1, inst)] = 1.0;
        }
        else
        {
            inst->solution[xpos(i, preds[i], j, inst)] = 1.0; 
        }
    }

    // evaluate C and total stream
    for (i = 0; i < inst->n_turbines; ++i)
    {
        if (inst->solution[fpos(i, inst->subs_index, inst)] > EPSILON_MED)
        {
            ++in_deg_cnt;
        }
    }

    if (in_deg_cnt > inst->subs_cables)
    {
        score += 10e12;
    }

    // evaluate streams
    for (i = 0; i < inst->n_turbines; ++i)
    {
        for (j = 0; j < inst->n_turbines; ++j)
        {
            if (inst->solution[fpos(i, j, inst)] > EPSILON_SMALL)
            {
                if ((inst->solution[fpos(i, j, inst)] - inst->cable_powers[inst->n_cables - 1]) > EPSILON_SMALL)
                {
                    score += 10e9 * (inst->solution[fpos(i, j, inst)] - inst->cable_powers[inst->n_cables - 1]);
                    ++stream_exceed_cnt;
                }
                else
                {
                    for (k = 0; k < inst->n_cables; ++k)
                    {
                        if ((inst->cable_powers[k] + EPSILON_SMALL) >= inst->solution[fpos(i, j, inst)])
                        {
                            score += inst->cable_costs[k] * hypot(inst->x_turb_coords[i] - inst->x_turb_coords[j], inst->y_turb_coords[i] - inst->y_turb_coords[j]);
                            break;
                        }
                    }
                }
            }

            if ((score - best_score) > EPSILON_MED)
            {
                return score;
            }
        }
    }

    if (inst->use_cross_table)
    {
       // evaluate cross with lookup table
        for (i = 0; i < inst->n_turbines; ++i)
        {
            for (j = 0; j < inst->n_turbines; ++j)
            {
                if ((inst->solution[ypos(i, j, inst)] + inst->solution[ypos(j, i, inst)]) < EPSILON_SMALL) continue;
                for (k = i+1; k < inst->n_turbines; ++k)
                {
                    for (h = 0; h < inst->n_turbines; ++h)
                    {
                        // true if (i,j) exists and (k,h) exists and (i,j) cross (k,h)
                        if ((inst->solution[ypos(i, j, inst)] > EPSILON_MED) && (inst->solution[ypos(k, h, inst)] > EPSILON_MED) && get_cross_table(inst, inst->cross_table, i, j, k, h))
                        {
                            score += 10e9;
                            ++cross_cnt;
                        }
                        if ((score - best_score) > EPSILON_MED)
                        {
                            return score;
                        }
                    }
                }
            }
        } 
    } else
    {
        // evaluate cross with lookup table
        for (i = 0; i < inst->n_turbines; ++i)
        {
            for (j = 0; j < inst->n_turbines; ++j)
            {
                if ((inst->solution[ypos(i, j, inst)] + inst->solution[ypos(j, i, inst)]) < EPSILON_SMALL) continue;
                for (k = i+1; k < inst->n_turbines; ++k)
                {
                    for (h = 0; h < inst->n_turbines; ++h)
                    {
                        // true if (i,j) exists and (k,h) exists and (i,j) cross (k,h)
                        if ((inst->solution[ypos(i, j, inst)] > EPSILON_MED) && (inst->solution[ypos(k, h, inst)] > EPSILON_MED) && is_crossing(inst->x_turb_coords[i], inst->y_turb_coords[i], inst->x_turb_coords[j], inst->y_turb_coords[j], inst->x_turb_coords[k], inst->y_turb_coords[k], inst->x_turb_coords[h], inst->y_turb_coords[h]))
                        {
                            score += 10e9;
                            ++cross_cnt;
                        }
                        if ((score - best_score) > EPSILON_MED)
                        {
                            return score;
                        }
                    }
                }
            }
        }
    }

    // evaluate cross with lookup table
    /*for (i = 0; i < inst->n_turbines; ++i)
    {
        for (j = 0; j < inst->n_turbines; ++j)
        {
            if ((inst->solution[ypos(i, j, inst)] + inst->solution[ypos(j, i, inst)]) < EPSILON_SMALL) continue;

            for (k = i+1; k < inst->n_turbines; ++k)
            {
                for (h = 0; h < inst->n_turbines; ++h)
                {
                    // true if (i,j) exists and (k,h) exists and (i,j) cross (k,h)
                    if ((inst->solution[ypos(i, j, inst)] > EPSILON_MED) && (inst->solution[ypos(k, h, inst)] > EPSILON_MED) && get_cross_table(inst, inst->cross_table, i, j, k, h))
                    {
                        score += 10e9;
                        ++cross_cnt;
                    }

                    if ((score - best_score) > EPSILON_MED)
                    {
                        return score;
                    }
                }
            }
        }
    }*/
    
    if (log)
    {
        printf("            in_deg =            %d\n", in_deg_cnt);
        printf("            tot_stream =        %lf\n", inst->streams[0]);
        printf("            stream_exceed_cnt = %d\n", stream_exceed_cnt);
        printf("            cross_cnt =         %d\n", cross_cnt);
    }

    return score;
}

// evaluate curr preds with prev preds, changing edge (a,b) to (a,c)
/*double evaluate_preds_param(instance *inst, int *preds, int a, int b, int c, double best_score, int log)
{
    int cnt, edges;
    int i, j, k, h, in_deg_cnt = 0, cross_cnt = 0, stream_exceed_cnt = 0;
    double score = 0.0;
    
    // clear and initialize structures
    // streams[i] is the input stream in the i-th turbine
    for (i = 0; i < inst->ncols; ++i) inst->solution[i] = 0.0;

    for (i = 0; i < inst->n_turbines; ++i)
    {
        inst->streams[i] = 0.0;
        inst->preds_copy[i] = preds[i];
        inst->flags[i] = 1;
    }

    

    if ((inst->streams[0] - inst->turb_tot_power) > EPSILON_SMALL || (inst->streams[0] - inst->turb_tot_power) < -EPSILON_SMALL)
    {
        // not acceptable
        return DBL_MAX;
    }

    for (i = 0; i < inst->n_turbines; ++i)
    {
        if (i == preds[i]) continue;

        inst->solution[fpos(i, preds[i], inst)] = inst->streams[i] + 1;
        inst->solution[ypos(i, preds[i], inst)] = 1.0;

        for (j = 0; j < inst->n_cables; ++j)
        {
            if (inst->cable_powers[j] >= inst->streams[i] + 1) break;
        }

        if (inst->cable_powers[inst->n_cables - 1] < inst->streams[i] + 1)
        {
            inst->solution[xpos(i, preds[i], inst->n_cables - 1, inst)] = 1.0;
        }
        else
        {
            inst->solution[xpos(i, preds[i], j, inst)] = 1.0; 
        }
    }

    // evaluate C and total stream
    for (i = 0; i < inst->n_turbines; ++i)
    {
        if (inst->solution[fpos(i, inst->subs_index, inst)] > EPSILON_MED)
        {
            ++in_deg_cnt;
        }
    }

    if (in_deg_cnt > inst->subs_cables)
    {
        score += 10e12;
    }

    // evaluate streams
    for (i = 0; i < inst->n_turbines; ++i)
    {
        for (j = 0; j < inst->n_turbines; ++j)
        {
            if (inst->solution[fpos(i, j, inst)] > EPSILON_SMALL)
            {
                if ((inst->solution[fpos(i, j, inst)] - inst->cable_powers[inst->n_cables - 1]) > EPSILON_SMALL)
                {
                    score += 10e9 * (inst->solution[fpos(i, j, inst)] - inst->cable_powers[inst->n_cables - 1]);
                    ++stream_exceed_cnt;
                }
                else
                {
                    for (k = 0; k < inst->n_cables; ++k)
                    {
                        if ((inst->cable_powers[k] + EPSILON_SMALL) >= inst->solution[fpos(i, j, inst)])
                        {
                            score += inst->cable_costs[k] * hypot(inst->x_turb_coords[i] - inst->x_turb_coords[j], inst->y_turb_coords[i] - inst->y_turb_coords[j]);
                            break;
                        }
                    }
                }
            }

            if ((score - best_score) > EPSILON_MED)
            {
                return score;
            }
        }
    }

    // evaluate cross
    for (i = 0; i < inst->n_turbines; ++i)
    {
        for (j = 0; j < inst->n_turbines; ++j)
        {
            if ((inst->solution[ypos(i, j, inst)] + inst->solution[ypos(j, i, inst)]) < EPSILON_SMALL) continue;

            for (k = i+1; k < inst->n_turbines; ++k)
            {
                for (h = 0; h < inst->n_turbines; ++h)
                {
                    // true if (i,j) exists and (k,h) exists and (i,j) cross (k,h)
                    if ((inst->solution[ypos(i, j, inst)] > EPSILON_MED) && (inst->solution[ypos(k, h, inst)] > EPSILON_MED) && is_crossing(inst->x_turb_coords[i], inst->y_turb_coords[i], inst->x_turb_coords[j], inst->y_turb_coords[j], inst->x_turb_coords[k], inst->y_turb_coords[k], inst->x_turb_coords[h], inst->y_turb_coords[h]))
                    {
                        score += 10e9;
                        ++cross_cnt;
                    }

                    if ((score - best_score) > EPSILON_MED)
                    {
                        return score;
                    }
                }
            }
        }
    }

    
    if (log)
    {
        printf("            in_deg =            %d\n", in_deg_cnt);
        printf("            tot_stream =        %lf\n", inst->streams[0]);
        printf("            stream_exceed_cnt = %d\n", stream_exceed_cnt);
        printf("            cross_cnt =         %d\n", cross_cnt);
    }

    return score;
}*/

// possibilimente da rimuovere
/*double evaluate_solution(instance *inst, double *solution, double best_score, int log)
{
    int i, j, k, h, in_deg_cnt = 0, cross_cnt = 0, stream_exceed_cnt = 0;
    double tot_stream = 0.0;
    double score = 0.0;

    // evaluate C and total stream
    for (i = 0; i < inst->n_turbines; ++i)
    {
        if (solution[fpos(i, inst->subs_index, inst)] > EPSILON_MED)
        {
            ++in_deg_cnt;
            tot_stream += solution[fpos(i, inst->subs_index, inst)];
        }
    }

    if (in_deg_cnt > inst-> subs_cables)
    {
        score += 10e12;
    }

    // WARN: si suppone che l'energia prodotta da ciascuna turbina sia 1.0
    if ((tot_stream - (inst->n_turbines - 1)) > EPSILON_SMALL || (tot_stream - (inst->n_turbines - 1)) < -EPSILON_SMALL)
    {
        // not acceptable
        return DBL_MAX;
    }

    // evaluate streams
    for (i = 0; i < inst->n_turbines; ++i)
    {
        for (j = 0; j < inst->n_turbines; ++j)
        {
            if (inst->solution[fpos(i, j, inst)] > EPSILON_SMALL)
            {
                if ((inst->solution[fpos(i, j, inst)] - inst->cable_powers[inst->n_cables - 1]) > EPSILON_SMALL)
                {
                    score += 10e9 * (inst->solution[fpos(i, j, inst)] - inst->cable_powers[inst->n_cables - 1]);
                    ++stream_exceed_cnt;
                }
                else
                {
                    for (k = 0; k < inst->n_cables; ++k)
                    {
                        if ((inst->cable_powers[k] + EPSILON_SMALL) >= inst->solution[fpos(i, j, inst)])
                        {
                            score += inst->cable_costs[k] * hypot(inst->x_turb_coords[i] - inst->x_turb_coords[j], inst->y_turb_coords[i] - inst->y_turb_coords[j]);
                            break;
                        }
                    }
                }
            }

            if ((score - best_score) > EPSILON_MED)
            {
                return score;
            }
        }
    }

    if (inst->use_cross_table)
    {
       // evaluate cross with lookup table
        for (i = 0; i < inst->n_turbines; ++i)
        {
            for (j = 0; j < inst->n_turbines; ++j)
            {
                if ((inst->solution[ypos(i, j, inst)] + inst->solution[ypos(j, i, inst)]) < EPSILON_SMALL) continue;
                for (k = i+1; k < inst->n_turbines; ++k)
                {
                    for (h = 0; h < inst->n_turbines; ++h)
                    {
                        // true if (i,j) exists and (k,h) exists and (i,j) cross (k,h)
                        if ((inst->solution[ypos(i, j, inst)] > EPSILON_MED) && (inst->solution[ypos(k, h, inst)] > EPSILON_MED) && get_cross_table(inst, inst->cross_table, i, j, k, h))
                        {
                            score += 10e9;
                            ++cross_cnt;
                        }
                        if ((score - best_score) > EPSILON_MED)
                        {
                            return score;
                        }
                    }
                }
            }
        } 
    } else
    {
        // evaluate cross with lookup table
        for (i = 0; i < inst->n_turbines; ++i)
        {
            for (j = 0; j < inst->n_turbines; ++j)
            {
                if ((inst->solution[ypos(i, j, inst)] + inst->solution[ypos(j, i, inst)]) < EPSILON_SMALL) continue;
                for (k = i+1; k < inst->n_turbines; ++k)
                {
                    for (h = 0; h < inst->n_turbines; ++h)
                    {
                        // true if (i,j) exists and (k,h) exists and (i,j) cross (k,h)
                        if ((inst->solution[ypos(i, j, inst)] > EPSILON_MED) && (inst->solution[ypos(k, h, inst)] > EPSILON_MED) && is_crossing(inst->x_turb_coords[i], inst->y_turb_coords[i], inst->x_turb_coords[j], inst->y_turb_coords[j], inst->x_turb_coords[k], inst->y_turb_coords[k], inst->x_turb_coords[h], inst->y_turb_coords[h]))
                        {
                            score += 10e9;
                            ++cross_cnt;
                        }
                        if ((score - best_score) > EPSILON_MED)
                        {
                            return score;
                        }
                    }
                }
            }
        }
    }
    
    if (log)
    {
        printf("            in_deg=             %d\n", in_deg_cnt);
        printf("            tot_stream=         %lf\n", tot_stream);
        printf("            stream_exceed_cnt = %d\n", stream_exceed_cnt);
        printf("            cross_cnt=          %d\n", cross_cnt);
    }

    return score;
}*/

int heur_VNS_launcher(instance *inst) {
    int j, test_cnt = 0;
    double score, best_score = DBL_MAX;

    inst->ncols = inst->n_turbines * inst->n_turbines * (inst->n_cables + 2);
    
    // strutture locali
    int *preds = (int *) calloc(sizeof(int), inst->n_turbines);
    int *best_preds = (int *) calloc(sizeof(int), inst->n_turbines);
    
    // strutture per le funzioni
    inst->flags = (int *) calloc(sizeof(int), inst->n_turbines);
    inst->preds_copy = (int *) calloc(sizeof(int), inst->n_turbines);
    
    inst->solution = (double *) calloc(sizeof(double), inst->ncols);
    inst->streams = (double *) calloc(sizeof(double), inst->n_turbines);
    inst->xstar = (double *) calloc(sizeof(double), inst->ncols);
    
    inst->fstart = 0;
    inst->xstart = inst->n_turbines * inst->n_turbines;
    inst->ystart = inst->n_turbines * inst->n_turbines * (inst->n_cables + 1);

    inst->cross_table = generate_cross_table(inst);

    do
    {
        //printf("%s Launching test %d\n", get_log(), test_cnt);
        srand(test_cnt + inst->randomseed);
        revisited_dijkstra(inst, preds);
        
        score = VNS(inst, preds, DBL_MAX);

        if (score < best_score)
        {
            for (j = 0; j < inst->n_turbines; ++j) best_preds[j] = preds[j]; 

            best_score = score;
            printf("%s time_elapsed=%lf Solution found with score=%.4e\n", get_log(), get_time_elapsed(inst), score);
        }

        test_cnt++;
    } while (!is_time_limit_expired(inst));

    score = evaluate_preds(inst, best_preds, DBL_MAX, 1);
    printf("\n%s Best obj-val=%lf\n", get_log(), score);
    printf("%s Best integer=%lf\n", get_log(), score);

    printf("%s time_elapsed=%lf final_score=%.4e\n", get_log(), get_time_elapsed(inst), score);

    preds_to_solution(inst, best_preds, inst->xstar);
    
    free(preds);
    free(best_preds);

    free(inst->flags);
    free(inst->preds_copy);
    free(inst->solution);
    free(inst->streams);
    free(inst->cross_table);
    
    return 0;
}

int heur_simulated_annealing_launcher(instance *inst) {
    int j, test_cnt = 0;
    double score, best_score = DBL_MAX;

    inst->ncols = inst->n_turbines * inst->n_turbines * (inst->n_cables + 2);
    
    // strutture locali
    int *preds = (int *) calloc(sizeof(int), inst->n_turbines);
    int *best_preds = (int *) calloc(sizeof(int), inst->n_turbines);
    
    // strutture per le funzioni
    inst->flags = (int *) calloc(sizeof(int), inst->n_turbines);
    inst->preds_copy = (int *) calloc(sizeof(int), inst->n_turbines);
    
    inst->solution = (double *) calloc(sizeof(double), inst->ncols);
    inst->streams = (double *) calloc(sizeof(double), inst->n_turbines);
    inst->xstar = (double *) calloc(sizeof(double), inst->ncols);
    
    inst->fstart = 0;
    inst->xstart = inst->n_turbines * inst->n_turbines;
    inst->ystart = inst->n_turbines * inst->n_turbines * (inst->n_cables + 1);

    inst->cross_table = generate_cross_table(inst);

    do
    {
        printf("%s Launching test %d\n", get_log(), test_cnt);
        srand(test_cnt + inst->randomseed);
        //dijkstra(inst, preds);
        revisited_dijkstra(inst, preds);
        
        score = simulated_annealing(inst, preds, best_score);
        //printf("SCORE=%lf, BEST_SCORE=%lf\n", score, best_score);
        
        if (score < best_score)
        {
            //printf("ACCETTO\n");
            for (j = 0; j < inst->n_turbines; ++j) best_preds[j] = preds[j];
            best_score = score;
            //printf("%s time_elapsed=%lf Solution found with score=%.4e\n", get_log(), get_time_elapsed(inst), score);
        }
        printf("%s time_elapsed=%lf Solution found with score=%.4e\n", get_log(), get_time_elapsed(inst), score);
        
        test_cnt++;
    } while (!is_time_limit_expired(inst));

    //score = evaluate_preds(inst, best_preds, DBL_MAX, 1);
    printf("\n%s Best obj-val=%lf\n", get_log(), best_score);
    printf("%s Best integer=%lf\n", get_log(), best_score);

    printf("%s time_elapsed=%lf final_score=%.4e\n", get_log(), get_time_elapsed(inst), best_score);

    preds_to_solution(inst, best_preds, inst->xstar);

    free(preds);
    free(best_preds);

    free(inst->flags);
    free(inst->preds_copy);
    free(inst->solution);
    free(inst->streams);
    free(inst->cross_table);
    
    return 0;
}

int revisited_dijkstra(instance *inst, int *preds) {
    int i, j, h, k;
    int c = 10; // numero cluster che finiranno dentro la substation
    
    double *costs = (double *) calloc(sizeof(double), inst->n_turbines * inst->n_turbines);
    int *flags = (int *) calloc(sizeof(int), inst->n_turbines);
    int *temp = (int *) calloc(sizeof(int), inst->n_turbines);

    //Costruisco matrice di adiacenza

    for (i = 0; i < inst->n_turbines; ++i)
    {
        for (j = i + 1; j < inst->n_turbines; ++j)
        {
            double dist = hypot(inst->x_turb_coords[i] - inst->x_turb_coords[j], inst->y_turb_coords[i] - inst->y_turb_coords[j]);
            
            costs[i*inst->n_turbines + j] = dist;
            costs[j*inst->n_turbines + i] = dist;
        }
    }

    //inizializzo substation
    flags[inst->subs_index] = 1;
    preds[inst->subs_index] = inst->subs_index;

    h = 0;
    //Collego il massimo numero di turbine alla substation
    while(h < c)
    {
        k = rand() % (inst->n_turbines-1) +1; //una turbina a caso tranne la substation

        if(flags[k] == 0)
        {
            flags[k] = 1;
            preds[k] = inst->subs_index;
            temp[h] = k;
            h++;
        }   
    }

    //Metto come predecessori degli archi restanti gli archi appena scelti e presenti in temp
    for (i = 0; i < inst->n_turbines; ++i)
    {
        double min_dist = DBL_MAX;
        if(flags[i] == 0)
        {
            // Controllo la turbina più vicina tra quelle in temp
            for(j = 0; j < c; ++j)
            {
                if(min_dist > costs[i*inst->n_turbines + temp[j]])
                {
                    min_dist = costs[i*inst->n_turbines + temp[j]];
                    preds[i] = temp[j]; 
                    flags[i] = 1;
                }
            }
        }
    }

    free(flags);
    free(costs);
    free(temp);

    return 0;
}

double simulated_annealing(instance *inst, int *preds, double c) {
	double T = 1.0;                // Initial and final temperature 
    double Tmin = 0.00001;              // Temperature at which iteration terminates 
    double alpha = 0.9;             // Decrease in temperature 
    double val2;
    int numIterations = 10;       // Number of iterations of annealing before decreasing temperature 

    int i, j, w, x, y, z;
    double cost, min_cost = c;

    int *best_preds = (int *) calloc(inst->n_turbines, sizeof(int));
    for (j = 0; j < inst->n_turbines; ++j) {
        best_preds[j] = preds[j];
    }

    /* La nostra soluzione di partenza ci è data da preds */
    /* min_cost è il valore della migliore soluzione */

    // Continues annealing until reaching minimum 
    // temperature 
    while (T > Tmin)
    {
        for (i = 0; i < numIterations; ++i)
        { 
            w = (int) ((rand() / ((double) RAND_MAX)) * inst->n_turbines);
            x = (int) ((rand() / ((double) RAND_MAX)) * inst->n_turbines);
            y = (int) ((rand() / ((double) RAND_MAX)) * inst->n_turbines);
            z = (int) ((rand() / ((double) RAND_MAX)) * inst->n_turbines);
            preds[w] = x;
            preds[y] = z;

            cost = two_opt_move(inst, preds, min_cost*(1 + T));
            //cost = two_opt_move(inst, preds, DBL_MAX);
            //cost = evaluate_preds(inst, preds, DBL_MAX, 0);
            
            //val1 = pow(M_E, - (1 - min_cost / cost));
            //val2 = pow(M_E, - (rand() / ((double) RAND_MAX)) * T);
            //val1 = (1 - (min_cost / cost));
            val2 = (rand() / ((double) RAND_MAX));
            //printf("val2=%lf, T=%lf\n", val2, T);
            if ((cost < min_cost) || val2 < T) {
            //if (cost < min_cost) {
                // accept
                for (j = 0; j < inst->n_turbines; ++j) {
                    best_preds[j] = preds[j];
                }
                min_cost = cost;
            } else {
                // reject
                for (j = 0; j < inst->n_turbines; ++j) {
                    preds[j] = best_preds[j];
                }
            }

            if (is_time_limit_expired(inst)) break;
        } 
  
        T *= alpha; // Decreases T, cooling phase 

        //printf("%s     Solution found with min_cost = %.4e\n", get_log(), min_cost);
        printf("%s time_elapsed=%4.4lf Solution found with cost=%.4e\n", get_log(), get_time_elapsed(inst), min_cost);
        
        if (inst->time_limit_expired) break;
    }

    evaluate_preds(inst, best_preds, DBL_MAX, 0);
    free(best_preds);
    return min_cost;
}
