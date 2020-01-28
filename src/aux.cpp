
#include "wf.h"
#include "aux.h"


// Verifica se i segmenti che hanno per estremi i punti (x1, y1), (x2, y2) e
// (x3, y3), (x4, y4) si intersecano o no.
int is_crossing(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4) {
    double A, B, C, D, E, F, det, lambda, mu;
    A = x1 - x2;
    B = x4 - x3;
    C = y1 - y2;
    D = y4 - y3;

    E = x4 - x2;
    F = y4 - y2;

    det = A * D - B * C;
    if (fabs(det) < EPSILON_SMALL)
        return 0; //no crossing
    else
    {
        lambda = (D * E - B * F) / det;
        mu = (A * F - C * E) / det;

        if ((lambda > EPSILON_MED) && (lambda < (1.0 - EPSILON_MED)) &&
            (mu > EPSILON_MED) && (mu < (1.0 - EPSILON_MED)))

        //if ((lambda > EPSILON_SMALL) && (lambda < (1.0 - EPSILON_SMALL)) &&
        //    (mu > EPSILON_SMALL) && (mu < (1.0 - EPSILON_SMALL)))

        {
            return 1; //crossing
        }
        return 0; //no crossing
    }
}

int no_cross_separation(instance *inst, double *xstar, int i, int j, int k, int h, int thread, int table) {
    int cross;
    double cable1, cable2;

    if (table)
    {
        cross = get_cross_table(inst, inst->cross_table, i, j, k, h);
    }
    else
    {
        cross = is_crossing(inst->x_turb_coords[i], inst->y_turb_coords[i], inst->x_turb_coords[j], inst->y_turb_coords[j], inst->x_turb_coords[k], inst->y_turb_coords[k], inst->x_turb_coords[h], inst->y_turb_coords[h]);
    }

    switch (inst->model_type)
    {
    case 2: // Model solved with loop method
        cable1 = inst->xstar[y_pos(i, j, inst)];
        cable2 = inst->xstar[y_pos(k, h, inst)];

        return cable1 > EPSILON_MED && cable2 > EPSILON_MED && cross == 1;
        break;
    case 3: // Model solved using callbacks

        if (inst->multi)
        {
            cable1 = inst->xstar_vector[thread][y_pos(i, j, inst)];
            cable2 = inst->xstar_vector[thread][y_pos(k, h, inst)];
        }
        else
        {
            cable1 = xstar[y_pos(i, j, inst)];
            cable2 = xstar[y_pos(k, h, inst)];
        }

        return cable1 > EPSILON_MED && cable2 > EPSILON_MED && cross == 1;
        break;
    default:
        return cross;
        break;
    }
}

int compute_no_cross_cuts(instance *inst, CPXCENVptr env, CPXLPptr lp, int thread) {
    char less_equal = 'L';
    int h, i, j, k, ncols, lastrow;
    int izero = 0;
    int nzcnt = 0;
    int ncuts = 0;
    double one = 1.0;

    ncols = CPXgetnumcols(env, lp);
    int *rmatind = (int *)calloc(ncols, sizeof(int));
    double *rmatval = (double *)calloc(ncols, sizeof(double));

    for (i = 0; i < inst->n_turbines; ++i)
    {
        for (j = i + 1; j < inst->n_turbines; ++j)
        {
            for (k = 0; k < inst->n_turbines; ++k)
            {
                if (k == i || k == j)
                    continue;

                // LOOP TASK: non serve controllare: nella soluzione trovata da CPLEX gli archi (i,j) e
                // (j,i) non sono stati scelti.
                if (inst->model_type == LOOP_TASK_MODEL_TYPE && (inst->xstar[y_pos(i, j, inst)] + inst->xstar[y_pos(j, i, inst)] < EPSILON_SMALL))
                    continue;

                nzcnt = 0;
                for (h = 0; h < inst->n_turbines; ++h)
                {
                    if (h == i || h == j || h == k)
                        continue;

                    if (no_cross_separation(inst, NULL, i, j, k, h, thread, inst->use_cross_table))
                    {
                        rmatind[nzcnt] = y_pos(k, h, inst);
                        rmatval[nzcnt] = 1.0;
                        nzcnt++;
                    }
                }

                if (nzcnt > 0)
                {
                    rmatind[nzcnt] = y_pos(i, j, inst);
                    rmatval[nzcnt] = 1.0;
                    nzcnt++;

                    rmatind[nzcnt] = y_pos(j, i, inst);
                    rmatval[nzcnt] = 1.0;
                    nzcnt++;

                    if (inst->model_type == LAZY_CONST_MODEL_TYPE)
                    {
                        // LAZY CONSTRAINT
                        if (CPXaddlazyconstraints(env, lp, 1, nzcnt, &one, &less_equal, &izero, rmatind, rmatval, NULL))
                            print_error("Wrong CPXaddlazyconstraints() at no-cross constraints.");
                    }

                    if (inst->model_type == NEW_ROWS_MODEL_TYPE)
                    {
                        // NEW ROW
                        lastrow = CPXgetnumrows(env, lp);
                        if (CPXnewrows(env, lp, 1, &one, &less_equal, NULL, NULL))
                            print_error("Wrong CPXnewrows() at no-cross constraints.");
                        for (h = 0; h < nzcnt; h++)
                        {
                            if (CPXchgcoef(env, lp, lastrow, rmatind[h], rmatval[h]))
                                print_error("Wrong CPXchgcoef() at no-cross constraints.");
                        }
                    }

                    if (inst->model_type == LOOP_TASK_MODEL_TYPE)
                    {
                        lastrow = CPXgetnumrows(env, lp);
                        CPXnewrows(env, lp, 1, &one, &less_equal, NULL, NULL);
                        for (h = 0; h < nzcnt; h++)
                        {
                            CPXchgcoef(env, lp, lastrow, rmatind[h], rmatval[h]);
                        }
                    }

                    ncuts++;
                }
            }
        }
    }

    free(rmatind);
    free(rmatval);
    return ncuts;
}

int CPXPUBLIC callback_cross_multi(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p) {
    int i, j, k, h, nzcnt, mythread;
    int ncuts = 0;

    // get mythread number
    CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MY_THREAD_NUM, &mythread);

    *useraction_p = CPX_CALLBACK_DEFAULT;
    instance *inst = (instance *)cbhandle; // casting of cbhandle

    if (inst->debug)
    {
        printf("message from thread=%d\n", mythread);
    }

    if (CPXgetcallbacknodex(env, cbdata, wherefrom, inst->xstar_vector[mythread], 0, inst->ncols - 1))
        return 1; // y = current y from CPLEX-- y starts from position 0

    for (i = 0; i < inst->n_turbines; ++i)
    {
        for (j = i + 1; j < inst->n_turbines; ++j)
        {
            if (inst->xstar_vector[mythread][y_pos(i, j, inst)] + inst->xstar_vector[mythread][y_pos(j, i, inst)] < EPSILON_SMALL)
                continue;

            for (k = 0; k < inst->n_turbines; ++k)
            {
                nzcnt = 0;
                for (h = 0; h < inst->n_turbines; ++h)
                {
                    if (h == i || h == j || h == k)
                        continue;

                    if (no_cross_separation(inst, NULL, i, j, k, h, mythread, inst->use_cross_table))
                    {
                        inst->indeces_vector[mythread][nzcnt] = y_pos(k, h, inst);
                        inst->values_vector[mythread][nzcnt] = 1.0;

                        nzcnt++;
                    }
                }
                if (nzcnt > 0)
                {
                    inst->indeces_vector[mythread][nzcnt] = y_pos(i, j, inst);
                    inst->values_vector[mythread][nzcnt] = 1.0;
                    nzcnt++;

                    inst->indeces_vector[mythread][nzcnt] = y_pos(j, i, inst);
                    inst->values_vector[mythread][nzcnt] = 1.0;
                    nzcnt++;

                    if (CPXcutcallbackadd(env, cbdata, wherefrom, nzcnt, 1.0, 'L', inst->indeces_vector[mythread], inst->values_vector[mythread], 0))
                        print_error("USER_separation: CPXcutcallbackadd error");
                    ncuts++;
                }
            }
        }
    }

    if (ncuts >= 1)
        *useraction_p = CPX_CALLBACK_SET; // tell CPLEX that cuts have been created
    return 0;                             // return 1 would mean error --> abort Cplex's execution
}

int CPXPUBLIC callback_cross_single(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p) {
    int i, j, k, h, nzcnt, mythread;
    int ncuts = 0;

    // get mythread number
    CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MY_THREAD_NUM, &mythread);

    *useraction_p = CPX_CALLBACK_DEFAULT;
    instance *inst = (instance *)cbhandle; // casting of cbhandle

    if (inst->debug)
    {
        printf("message from thread=%d\n", mythread);
    }

    double *xstar = (double *)calloc(inst->ncols, sizeof(double));
    int *indeces = (int *)calloc(inst->ncols, sizeof(int));
    double *values = (double *)calloc(inst->ncols, sizeof(double));

    if (CPXgetcallbacknodex(env, cbdata, wherefrom, xstar, 0, inst->ncols - 1))
        return 1; // y = current y from CPLEX-- y starts from position 0

    for (i = 0; i < inst->n_turbines; ++i)
    {
        for (j = 0; j < inst->n_turbines; ++j)
        {
            if (xstar[y_pos(i, j, inst)] + xstar[y_pos(j, i, inst)] < EPSILON_SMALL)
                continue;

            for (k = 0; k < inst->n_turbines; ++k)
            {
                nzcnt = 0;
                for (h = 0; h < inst->n_turbines; ++h)
                {
                    //if ((xstar[y_pos(i, j, inst)] > EPSILON_MED) && (xstar[y_pos(k, h, inst)] > EPSILON_MED) && get_cross_table(inst, inst->cross_table, i, j, k, h))
                    if ((xstar[y_pos(i, j, inst)] > EPSILON_MED) && (xstar[y_pos(k, h, inst)] > EPSILON_MED) && no_cross_separation(inst, xstar, i, j, k, h, mythread, inst->use_cross_table))

                    {
                        indeces[nzcnt] = y_pos(k, h, inst);
                        values[nzcnt] = 1.0;
                        nzcnt++;
                    }
                }
                if (nzcnt > 0)
                {
                    indeces[nzcnt] = y_pos(i, j, inst);
                    values[nzcnt] = 1.0;
                    nzcnt++;

                    indeces[nzcnt] = y_pos(j, i, inst);
                    values[nzcnt] = 1.0;
                    nzcnt++;

                    if (CPXcutcallbackadd(env, cbdata, wherefrom, nzcnt, 1.0, 'L', indeces, values, 0))
                        print_error("USER_separation: CPXcutcallbackadd error");

                    ncuts++;
                }
            }
        }
    }

    free(xstar);
    free(indeces);
    free(values);
    if (ncuts >= 1)
        *useraction_p = CPX_CALLBACK_SET; // tell CPLEX that cuts have been created
    return 0;                             // return 1 would mean error --> abort Cplex's execution
}

char *generate_cross_table(instance *inst) {
    int i, j, k, h;

    // allocazione completa
    //int byte_num = ceil(pow(inst->n_turbines, 4) / 8.0);
    //char *cross_table = (char *) calloc(sizeof(char), byte_num);

    /*  Allocazione parziale
        
        In questo tipo di allocazione la cross_table viene allocata solo in modo parziale
        (con notevole risparmio di memoria) grazie alla proprietà di invarianza degli indici
        (i,j) e (k,h), infatti, al fine del calcolo dei cross l'arco (i,j) == (j,i).

        Gli elementi della tabella non sono più quindi n^4 ma (n^4 - 2*n^3 + n^2) / 4 (con n
        numero delle turbine) poiché il lato della tabella quadrata vale n*(n-1) / 2.

        Per calcolare l'indice di riga (o di colonna) dati i due indici delle turbine (i,j) vale
        la formula:

        off = i*(n-1) + j - i*(i+1) / 2 - 1

        con 0 <= i < j < n.

        Per calcolare l'offset della struttura lineare si usano le formule note per le matrici
        quadrate (ricordando il valore del lato).
    */
    int byte_num = ceil(((pow(inst->n_turbines, 4) - 2 * pow(inst->n_turbines, 3) + pow(inst->n_turbines, 2)) / 4.0) / 8.0);
    char *cross_table = (char *)calloc(sizeof(char), byte_num);

    for (i = 0; i < inst->n_turbines; ++i)
    {
        for (j = 0; j < inst->n_turbines; ++j)
        {
            for (k = 0; k < inst->n_turbines; ++k)
            {
                for (h = 0; h < inst->n_turbines; ++h)
                {
                    if (is_crossing(inst->x_turb_coords[i], inst->y_turb_coords[i], inst->x_turb_coords[j], inst->y_turb_coords[j], inst->x_turb_coords[k], inst->y_turb_coords[k], inst->x_turb_coords[h], inst->y_turb_coords[h]))
                    {
                        set_cross_table(inst, cross_table, 1, i, j, k, h);
                    }
                }
            }
        }
    }

    return cross_table;
}

void set_cross_table(instance *inst, char *cross_table, char value, int i, int j, int k, int h) {
    int tmp, row_offset, col_offset, offset, bit_offset, byte_offset;

    // allocazione completa
    //offset = ((i * inst->n_turbines + j) * inst->n_turbines + k) * inst->n_turbines + h;

    // allocazione parziale
    if (i > j)
    {
        tmp = i;
        i = j;
        j = tmp;
    }
    else
    {
        if (i == j)
            return;
    }

    if (k > h)
    {
        tmp = k;
        k = h;
        h = tmp;
    }
    else
    {
        if (k == h)
            return;
    }

    row_offset = i * (inst->n_turbines - 1) + j - ((i * i + i) >> 1) - 1;
    col_offset = k * (inst->n_turbines - 1) + h - ((k * k + k) >> 1) - 1;
    offset = row_offset * ((inst->n_turbines * inst->n_turbines - inst->n_turbines) >> 1) + col_offset;

    // parte comune
    byte_offset = offset >> 3;
    bit_offset = offset % 8;

    if (value)
    {
        if (cross_table[byte_offset] & (1 << bit_offset))
            return;
        else
            cross_table[byte_offset] = cross_table[byte_offset] + (1 << bit_offset);
    }
    else
    {
        if (cross_table[byte_offset] & (1 << bit_offset))
            cross_table[byte_offset] = cross_table[byte_offset] - (1 << bit_offset);
        else
            return;
    }
}

int get_cross_table(instance *inst, char *cross_table, int i, int j, int k, int h) {
    // allocazione completa
    //int var = ((i * inst->n_turbines + j) * inst->n_turbines + k) * inst->n_turbines + h;

    // allocazione parziale
    int var;
    if (i > j)
    {
        var = i;
        i = j;
        j = var;
    }
    else
    {
        if (i == j)
            return 0;
    }

    if (k > h)
    {
        var = k;
        k = h;
        h = var;
    }
    else
    {
        if (k == h)
            return 0;
    }

    var = (i * (inst->n_turbines - 1) + j - ((i * i + i) >> 1) - 1) * ((inst->n_turbines * inst->n_turbines - inst->n_turbines) >> 1) + (k * (inst->n_turbines - 1) + h - ((k * k + k) >> 1) - 1);

    return (cross_table[var >> 3] & (1 << (var % 8))) ? 1 : 0;
}

int hard_fix(instance *inst, CPXENVptr env, CPXLPptr lp) {
    int i, j;
    char *lu;
    double *bd;
    int *indices;
    lu = (char *) calloc(inst->n_turbines * inst->n_turbines, sizeof(char));
    bd = (double *) calloc(inst->n_turbines * inst->n_turbines, sizeof(double));
    indices = (int *) calloc(inst->n_turbines * inst->n_turbines, sizeof(int));

    for (i = 0; i < inst->n_turbines; ++i)
    {
        for (j = 0; j < inst->n_turbines; ++j)
        {
            indices[inst->n_turbines * i + j] = inst->ystart + inst->n_turbines * i + j;

            if ((rand() / ((double)RAND_MAX)) < inst->hard_fix_prob)
            {
                bd[inst->n_turbines * i + j] = inst->xstar[inst->ystart + inst->n_turbines * i + j];
                lu[inst->n_turbines * i + j] = 'B';
            }

            else
            {
                bd[inst->n_turbines * i + j] = 0;
                lu[inst->n_turbines * i + j] = 'L';
            }
        }
    }

    CPXchgbds(env, lp, inst->n_turbines * inst->n_turbines, indices, lu, bd);

    free(lu);
    free(bd);
    free(indices);

    return 0;
}

int unfix_variables(instance *inst, CPXENVptr env, CPXLPptr lp) {
    int i, j;
    char *lu;
    double *bd;
    int *indices;
    lu = (char *) calloc(inst->n_turbines * inst->n_turbines, sizeof(char));
    bd = (double *) calloc(inst->n_turbines * inst->n_turbines, sizeof(double));
    indices = (int *) calloc(inst->n_turbines * inst->n_turbines, sizeof(int));

    for (i = 0; i < inst->n_turbines; ++i)
    {
        for (j = 0; j < inst->n_turbines; ++j)
        {
            indices[inst->n_turbines * i + j] = inst->ystart + inst->n_turbines * i + j;

            bd[inst->n_turbines * i + j] = 0;
            lu[inst->n_turbines * i + j] = 'L';
        }
    }

    CPXchgbds(env, lp, inst->n_turbines * inst->n_turbines, indices, lu, bd);

    for (i = 0; i < inst->n_turbines; ++i)
    {
        for (j = 0; j < inst->n_turbines; ++j)
        {
            indices[inst->n_turbines * i + j] = inst->ystart + inst->n_turbines * i + j;

            bd[inst->n_turbines * i + j] = 1;
            lu[inst->n_turbines * i + j] = 'U';
        }
    }

    CPXchgbds(env, lp, inst->n_turbines * inst->n_turbines, indices, lu, bd);

    free(lu);
    free(bd);
    free(indices);

    return 0;
}
