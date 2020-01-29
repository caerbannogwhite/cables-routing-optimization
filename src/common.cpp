
#include "common.hpp"

int comm_free_instance(instance *inst)
{
    free(inst->x_turb_coords);
    free(inst->y_turb_coords);
    free(inst->turb_powers);
    free(inst->cable_powers);
    free(inst->cable_costs);

    return 0;
}

int comm_read_input(instance *inst)
{
    char line[100];
    int index, line_counter;
    int x, y, p;
    double c, l;

    FILE *fin;

    // TURBINES
    fin = fopen(&inst->turbine_input_file[0], "r");
    if (fin == NULL)
        comm_print_error("Turbine input file not found.");

    line_counter = 0;

    // get number of turbine entries
    while (fgets(line, sizeof(line), fin) != NULL)
    {
        if (feof(fin))
            break;
        line_counter++;
    }
    inst->n_turbines = line_counter;
    fclose(fin);

    inst->x_turb_coords = (double *)calloc(inst->n_turbines, sizeof(double));
    inst->y_turb_coords = (double *)calloc(inst->n_turbines, sizeof(double));
    inst->turb_powers = (double *)calloc(inst->n_turbines, sizeof(double));

    // read the file and insert data
    fin = fopen(&inst->turbine_input_file[0], "r");
    for (index = 0; index < inst->n_turbines; index++)
    {
        if (!fscanf(fin, "%d %d %d", &x, &y, &p))
            comm_print_error("Wrong turbine file format.");

        inst->x_turb_coords[index] = (double)x;
        inst->y_turb_coords[index] = (double)y;
        inst->turb_powers[index] = (double)p;

        // looking for the substation, where power is -1
        if (p < -0.5)
        {
            inst->subs_index = index;
            inst->x_subs_coord = (double)x;
            inst->y_subs_coord = (double)y;
        }
        else
        {
            inst->turb_tot_power += p;
        }

        if (VERBOSE >= 1000)
            printf("x=%lf, y=%lf, power=%lf\n", inst->x_turb_coords[index], inst->y_turb_coords[index], inst->turb_powers[index]);
    }
    fclose(fin);

    // CABLES
    fin = fopen(&inst->cable_input_file[0], "r");
    if (fin == NULL)
        comm_print_error("Cable input file not found.");

    line_counter = 0;

    // get number of cable entries
    while (fgets(line, sizeof(line), fin) != NULL)
    {
        if (feof(fin))
            break;
        line_counter++;
    }
    inst->n_cables = line_counter;
    fclose(fin);

    inst->cable_powers = (double *)calloc(inst->n_cables, sizeof(double));
    inst->cable_costs = (double *)calloc(inst->n_cables, sizeof(double));

    // read the file and insert data
    fin = fopen(&inst->cable_input_file[0], "r");
    for (index = 0; index < inst->n_cables; index++)
    {
        if (!fscanf(fin, "%d %lf %lf", &p, &c, &l))
            comm_print_error("Wrong cable file format.");

        inst->cable_powers[index] = p;
        inst->cable_costs[index] = c;

        if (VERBOSE >= 1000)
            printf("power=%lf, cost=%lf\n", inst->cable_powers[index], inst->cable_costs[index]);
    }
    fclose(fin);

    return 0;
}

// Verifica se i segmenti che hanno per estremi i punti (x1, y1), (x2, y2) e
// (x3, y3), (x4, y4) si intersecano o no.
int comm_is_crossing(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
{
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

int comm_no_cross_separation(instance *inst, double *xstar, int i, int j, int k, int h, int thread, int table)
{
    int cross;
    double cable1, cable2;

    if (table)
    {
        cross = comm_get_cross_table(inst, inst->cross_table, i, j, k, h);
    }
    else
    {
        cross = comm_is_crossing(inst->x_turb_coords[i], inst->y_turb_coords[i], inst->x_turb_coords[j], inst->y_turb_coords[j], inst->x_turb_coords[k], inst->y_turb_coords[k], inst->x_turb_coords[h], inst->y_turb_coords[h]);
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

char *comm_generate_cross_table(instance *inst)
{
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
                    if (comm_is_crossing(inst->x_turb_coords[i], inst->y_turb_coords[i], inst->x_turb_coords[j], inst->y_turb_coords[j], inst->x_turb_coords[k], inst->y_turb_coords[k], inst->x_turb_coords[h], inst->y_turb_coords[h]))
                    {
                        comm_set_cross_table(inst, cross_table, 1, i, j, k, h);
                    }
                }
            }
        }
    }

    return cross_table;
}

void comm_set_cross_table(instance *inst, char *cross_table, char value, int i, int j, int k, int h)
{
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

int comm_get_cross_table(instance *inst, char *cross_table, int i, int j, int k, int h)
{
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

void comm_print_error(const char *err)
{
    time_t curtime;
    struct tm *tmp;

    curtime = time(NULL);
    tmp = gmtime(&curtime);

    fprintf(stderr, "\n[ERR %2d:%2d:%2d] %s\n", tmp->tm_hour, tmp->tm_min, tmp->tm_sec, err);
    fflush(NULL);
    exit(1);
}

int comm_is_time_limit_expired(instance *inst)
{
    clock_gettime(CLOCK_MONOTONIC, inst->time_end);

    double tspan = ((double)inst->time_end->tv_sec + 1.0e-9 * inst->time_end->tv_nsec) - ((double)inst->time_start->tv_sec + 1.0e-9 * inst->time_start->tv_nsec);
    if (tspan > inst->time_limit)
    {
        if (VERBOSE >= 100)
            printf("Time limit of %lf expired.\n", inst->time_limit);
        inst->time_limit_expired = 1;
        return 1;
    }
    return 0;
}

double comm_get_time_elapsed(instance *inst)
{
    clock_gettime(CLOCK_MONOTONIC, inst->time_end);
    double tspan = ((double)inst->time_end->tv_sec + 1.0e-9 * inst->time_end->tv_nsec) - ((double)inst->time_start->tv_sec + 1.0e-9 * inst->time_start->tv_nsec);
    return tspan;
}