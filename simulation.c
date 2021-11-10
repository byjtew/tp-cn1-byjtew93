#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define real float

// Solar constant W / m ^ 2
#define S0 1370

// Stefan - Boltzmann constant W / m2 / K4
#define SIGMA 0.00000005670367

// Temperature inertia (in years)
#define THETA 100.0

// Albedo, in this simulation albedo is considered constant
// in reality albedo can become lower with increased temperatures
// due to ice cap melting
#define ALBEDO 0.33

// Initial values

// Simulation starts in 2007
static const real t0 = 2007.0;
// Temperaturature in 2007 in K
static const real T0 = 288.45;
// CO2 concentration in 2007 in ppm
static const real CO20 = 370.0;

// User-defined:
static const unsigned minSteps = 16;

// Greenhouse gaz fraction
static real G(real T, real co2) {
    return 3.35071874e-03 * T + 3.20986702e-05 * co2 - 0.561593690144655;
}

static const real P_in = (1 - ALBEDO) * (S0 / 4.0);

static real P_out(real t, real T) {
    return T * T * T * T * SIGMA * (1.0 - G(T, t));
}

static real F(real t, real T) { return (P_in - P_out(t, T)) / 100.0; }

static real euler(real t_final, int steps, real* results) {
    real h = (t_final - t0) / steps;
    real temp = T0;
    real time = t0;
    real cur_CO2 = CO20;
    for (size_t i = 0; i < 2 * steps; i += 2) {
        temp += h * F(cur_CO2, temp);
        time += h;
        cur_CO2 *= 1.0000;
        results[i] = time;
        results[i + 1] = temp;
        //  printf("%f %f\n", time, temp);
        //  fprintf(fp, "%f %f\n", time, temp);
    }
    return temp;
}

static void errors(real* results, int steps, double result) {
    FILE *fpOut;
    real diff = 0.0;
    fpOut = fopen("errors.dat", "wb");
    for (size_t i = minSteps; i <= steps; i*=2) {
        diff = results[i-minSteps] > result ? results[i-minSteps] - result
                                     : result - results[i-minSteps];
        //printf("<> errors(): diff for %lu steps: %lf - %lf = %lf\n", i, result, results[i-minSteps], diff);
        fprintf(fpOut, "%lu %lf\n", i, diff);
    }
    fclose(fpOut);
}

int main(int argc, char** argv) {
    if (argc > 3 || argc < 2 || atoi(argv[1]) <= 10) {
        fprintf(
            stderr,
            "Usage: %s [steps > 10] [OPTIONAL: ANYTHING FOR ERROR TESTING]\n",
            argv[0]);
        exit(-1);
    }
    real* results;
    FILE* fp;
    size_t steps = atoi(argv[1]);

    if (argc > 2) {  // errors mesurement only
        printf(
            "Simple error measurement of the climate simulation with steps in "
            "range(%u, %lu)\n",
            minSteps, steps);
        real* tmp = malloc(sizeof(real) * steps * 2);
        results = malloc(sizeof(real) * steps);
        for (size_t i = minSteps; i <= steps; i*=2) {
            results[i - minSteps] = euler(2107, i, tmp);
        }
        errors(results, steps, 288.695113);
        free(tmp);
        printf("max steps = %lu; referenced result = %lf\nVisualize with: gnuplot plot-simulation_errors.gp\n", steps, 288.695113);

    } else if (argc == 2) {
        fp = fopen("reference.dat", "wb");
        printf("Simple climate simulation\n");
        results = malloc(sizeof(real) * steps * 2);
        real result = euler(2107, steps, results);
        for (size_t i = 0; i < steps * 2; i += 2) {
            fprintf(fp, "%lf %f\n", results[i], results[i + 1]);
            //printf("%lu %lf %f\n", i, results[i], results[i + 1]);
        }
        fclose(fp);
        printf("steps = %lu; final_result = %f\ngnuplot plot-simulation.gp\n", steps, result);
    }
    free(results);

    return 0;
}