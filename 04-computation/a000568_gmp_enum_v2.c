/*
 * a000568_gmp_enum_v2.c — A000568 via partition enumeration with GMP.
 *
 * Optimizations over v1:
 *   1. Incremental LCD/z computation: pass lcd_over_z as parameter,
 *      multiply by precomputed LCD/(k^m * m!) factor at each step.
 *      Avoids recomputing z_val from scratch and expensive division.
 *   2. Precomputed quotient table: qt[k][m] = LCD / (k^m * m!)
 *      for all valid (k,m) pairs. The incremental update is:
 *        new_lcd_over_z = lcd_over_z * qt[k][m] / LCD
 *      But that's still a division. Better approach:
 *      We carry z_inv_factor = product of (LCD / (k_i^m_i * m_i!))
 *      Actually simplest: carry z_val incrementally, divide once at leaf.
 *   3. Actually, the key insight: carry lcd_over_z incrementally.
 *      When adding part (k,m): lcd_over_z /= (k^m * m!)
 *      When removing: lcd_over_z *= (k^m * m!)
 *      This replaces leaf-level z computation + division with
 *      incremental multiply/divide at each recursion step.
 *
 * Compile:
 *   gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *       -o a000568_gmp_enum_v2 a000568_gmp_enum_v2.c -lgmp
 *
 * Author: opus-2026-03-08-S48
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>

#define MAX_N 1000
#define MAX_PARTS 500

/* GCD table */
static int gcd_tab[MAX_N + 1][MAX_N + 1];

static int gcd_int(int a, int b) {
    while (b) { int t = b; b = a % b; a = t; }
    return a;
}

static void init_gcd(int n) {
    for (int i = 1; i <= n; i += 2)
        for (int j = i; j <= n; j += 2) {
            int g = gcd_int(i, j);
            gcd_tab[i][j] = gcd_tab[j][i] = g;
        }
}

/* Odd parts list */
static int odd_parts[MAX_PARTS];
static int num_odd;

/* Precomputed factorials */
static mpz_t fact[MAX_N + 1];

/* LCD */
static mpz_t LCD;

/* Precomputed k^m * m! for all valid (k,m) pairs */
/* km_factor[k][m] = k^m * m! */
#define MAX_M_PER_K 502
static mpz_t km_factor[MAX_N + 1][MAX_M_PER_K];
static int km_factor_max_m[MAX_N + 1]; /* max m for each k */
static int km_inited[MAX_N + 1]; /* flag */

/* Accumulator */
static mpz_t total_sum;

/* Temporaries */
static mpz_t tmp_pow2, tmp_contrib;

/* Stack */
typedef struct { int k; int m; } PartEntry;
static PartEntry stack[MAX_PARTS];
static int stack_size = 0;

static long long part_count = 0;

/* lcd_over_z is passed by pointer for incremental update */
static void enumerate(int remaining, int max_part_idx,
                       long long t_val, mpz_t lcd_over_z)
{
    if (remaining == 0) {
        part_count++;
        /* contrib = lcd_over_z * 2^t_val */
        mpz_mul_2exp(tmp_contrib, lcd_over_z, (unsigned long)t_val);
        mpz_add(total_sum, total_sum, tmp_contrib);
        return;
    }

    mpz_t new_lcd_over_z;
    mpz_init(new_lcd_over_z);

    for (int pi = max_part_idx; pi >= 0; pi--) {
        int k = odd_parts[pi];
        if (k > remaining) continue;
        int max_m = remaining / k;

        for (int m = 1; m <= max_m; m++) {
            /* t update */
            long long dt_self = (long long)m * (m - 1) * k / 2 + (long long)m * (k - 1) / 2;
            long long dt_cross = 0;
            for (int s = 0; s < stack_size; s++) {
                dt_cross += (long long)m * stack[s].m * gcd_tab[k][stack[s].k];
            }
            long long new_t = t_val + dt_self + dt_cross;

            /* lcd_over_z update: divide by k^m * m! */
            mpz_tdiv_q(new_lcd_over_z, lcd_over_z, km_factor[k][m]);

            stack[stack_size].k = k;
            stack[stack_size].m = m;
            stack_size++;

            enumerate(remaining - m * k, pi - 1, new_t, new_lcd_over_z);

            stack_size--;
        }
    }

    mpz_clear(new_lcd_over_z);
}

void a000568(int n) {
    struct timespec t_start, t_end;
    clock_gettime(CLOCK_MONOTONIC, &t_start);

    if (n <= 1) {
        printf("a(%d) = 1\n", n);
        return;
    }

    init_gcd(n);

    num_odd = 0;
    for (int k = 1; k <= n; k += 2)
        odd_parts[num_odd++] = k;

    /* Init factorials */
    for (int i = 0; i <= n; i++)
        mpz_init(fact[i]);
    mpz_set_ui(fact[0], 1);
    for (int i = 1; i <= n; i++)
        mpz_mul_ui(fact[i], fact[i - 1], i);

    /* Precompute km_factor[k][m] = k^m * m! for each odd k */
    memset(km_inited, 0, sizeof(km_inited));
    for (int ki = 0; ki < num_odd; ki++) {
        int k = odd_parts[ki];
        int max_m = n / k;
        km_factor_max_m[k] = max_m;
        km_inited[k] = 1;
        for (int m = 1; m <= max_m; m++) {
            mpz_init(km_factor[k][m]);
            mpz_ui_pow_ui(km_factor[k][m], k, m);
            mpz_mul(km_factor[k][m], km_factor[k][m], fact[m]);
        }
    }

    /* Compute LCD */
    mpz_init_set_ui(LCD, 1);
    for (int k = 1; k <= n; k += 2) {
        int mk = n / k;
        mpz_mul(LCD, LCD, km_factor[k][mk]);
    }

    /* Init accumulators */
    mpz_init_set_ui(total_sum, 0);
    mpz_init(tmp_pow2);
    mpz_init(tmp_contrib);

    /* Enumerate starting with lcd_over_z = LCD (z=1 initially) */
    stack_size = 0;
    part_count = 0;
    enumerate(n, num_odd - 1, 0, LCD);

    /* Divide by LCD */
    mpz_t result;
    mpz_init(result);
    mpz_tdiv_q(result, total_sum, LCD);

    /* Verify */
    mpz_t remainder;
    mpz_init(remainder);
    mpz_tdiv_r(remainder, total_sum, LCD);
    if (mpz_sgn(remainder) != 0) {
        fprintf(stderr, "ERROR: not exactly divisible!\n");
    }
    mpz_clear(remainder);

    clock_gettime(CLOCK_MONOTONIC, &t_end);
    double elapsed = (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_nsec - t_start.tv_nsec) * 1e-9;

    size_t ndigits = mpz_sizeinbase(result, 10);
    printf("a(%d) = ", n);
    mpz_out_str(stdout, 10, result);
    printf("\n");
    printf("(%zu digits, %.1fs, %lld partitions)\n", ndigits, elapsed, part_count);

    /* Cleanup */
    mpz_clear(result);
    mpz_clear(total_sum);
    mpz_clear(LCD);
    mpz_clear(tmp_pow2);
    mpz_clear(tmp_contrib);
    for (int ki = 0; ki < num_odd; ki++) {
        int k = odd_parts[ki];
        for (int m = 1; m <= km_factor_max_m[k]; m++)
            mpz_clear(km_factor[k][m]);
    }
    for (int i = 0; i <= n; i++)
        mpz_clear(fact[i]);
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <n> [n2 n3 ...]\n", argv[0]);
        return 1;
    }

    for (int i = 1; i < argc; i++) {
        int n = atoi(argv[i]);
        a000568(n);
        fflush(stdout);
    }

    return 0;
}
