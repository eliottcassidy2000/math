/*
 * a000568_gmp_enum.c — A000568 via partition enumeration with GMP exact arithmetic.
 *
 * Replaces Python a000568_enum.py with C + GMP for 20-50x speedup.
 * Uses LCD-scaled integer accumulation (same algorithm as Python version).
 *
 * Compile:
 *   gcc -O3 -o a000568_gmp_enum a000568_gmp_enum.c -lgmp
 *
 * Usage:
 *   ./a000568_gmp_enum <n>
 *
 * Author: opus-2026-03-08-S47
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>

#define MAX_N 1000
#define MAX_PARTS 500  /* max distinct odd parts */

/* GCD table for odd numbers up to MAX_N */
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

/* Precomputed factorials as GMP integers */
static mpz_t fact[MAX_N + 1];

/* LCD = product of all possible z_lambda denominators */
static mpz_t LCD;

/* z-factor cache: z_factor[k][m] = k^m * m! */
/* We'll compute these on the fly to save memory */

/* Accumulator for the sum */
static mpz_t total_sum;

/* Temporary GMP variables (avoid repeated alloc) */
static mpz_t tmp_z_val, tmp_contrib, tmp_lcd_over_z, tmp_pow2, tmp_factor;

/* Partition stack */
typedef struct { int k; int m; } PartEntry;
static PartEntry stack[MAX_PARTS];
static int stack_size = 0;

/* Count partitions processed */
static long long part_count = 0;

static void enumerate(int n, int remaining, int max_part_idx,
                       long long t_val, /* current t value (fits in int64 for n < 500) */
                       int z_changed)    /* flag: z_val needs update */
{
    if (remaining == 0) {
        part_count++;

        /* Compute z_val = product of k^m * m! for each (k,m) in stack */
        mpz_set_ui(tmp_z_val, 1);
        for (int i = 0; i < stack_size; i++) {
            int k = stack[i].k;
            int m = stack[i].m;
            /* tmp_z_val *= k^m * m! */
            mpz_t km;
            mpz_init(km);
            mpz_ui_pow_ui(km, k, m);
            mpz_mul(tmp_z_val, tmp_z_val, km);
            mpz_mul(tmp_z_val, tmp_z_val, fact[m]);
            mpz_clear(km);
        }

        /* contrib = LCD / z_val * 2^t_val */
        mpz_tdiv_q(tmp_lcd_over_z, LCD, tmp_z_val);
        /* 2^t_val */
        mpz_set_ui(tmp_pow2, 1);
        mpz_mul_2exp(tmp_pow2, tmp_pow2, (unsigned long)t_val);
        /* contrib = lcd_over_z * pow2 */
        mpz_mul(tmp_contrib, tmp_lcd_over_z, tmp_pow2);
        /* accumulate */
        mpz_add(total_sum, total_sum, tmp_contrib);
        return;
    }

    for (int pi = max_part_idx; pi >= 0; pi--) {
        int k = odd_parts[pi];
        if (k > remaining) continue;
        int max_m = remaining / k;

        for (int m = 1; m <= max_m; m++) {
            /* Compute incremental t update */
            /* Self: m(m-1)k/2 + m(k-1)/2 */
            long long dt_self = (long long)m * (m - 1) * k / 2 + (long long)m * (k - 1) / 2;

            /* Cross: m * sum_{prev} m_prev * gcd(k, k_prev) */
            long long dt_cross = 0;
            for (int s = 0; s < stack_size; s++) {
                dt_cross += (long long)m * stack[s].m * gcd_tab[k][stack[s].k];
            }

            long long new_t = t_val + dt_self + dt_cross;

            stack[stack_size].k = k;
            stack[stack_size].m = m;
            stack_size++;

            enumerate(n, remaining - m * k, pi - 1, new_t, 1);

            stack_size--;
        }
    }
}

void a000568(int n) {
    struct timespec t_start, t_end;
    clock_gettime(CLOCK_MONOTONIC, &t_start);

    if (n <= 1) {
        printf("a(%d) = 1\n", n);
        return;
    }

    /* Init GCD table */
    init_gcd(n);

    /* Init odd parts */
    num_odd = 0;
    for (int k = 1; k <= n; k += 2)
        odd_parts[num_odd++] = k;

    /* Init factorials */
    for (int i = 0; i <= n; i++)
        mpz_init(fact[i]);
    mpz_set_ui(fact[0], 1);
    for (int i = 1; i <= n; i++)
        mpz_mul_ui(fact[i], fact[i - 1], i);

    /* Compute LCD = product_{k odd, k<=n} k^{floor(n/k)} * floor(n/k)! */
    mpz_init_set_ui(LCD, 1);
    for (int k = 1; k <= n; k += 2) {
        int mk = n / k;
        mpz_t kpow;
        mpz_init(kpow);
        mpz_ui_pow_ui(kpow, k, mk);
        mpz_mul(LCD, LCD, kpow);
        mpz_mul(LCD, LCD, fact[mk]);
        mpz_clear(kpow);
    }

    /* Init accumulators */
    mpz_init_set_ui(total_sum, 0);
    mpz_init(tmp_z_val);
    mpz_init(tmp_contrib);
    mpz_init(tmp_lcd_over_z);
    mpz_init(tmp_pow2);
    mpz_init(tmp_factor);

    /* Enumerate */
    stack_size = 0;
    part_count = 0;
    enumerate(n, n, num_odd - 1, 0, 0);

    /* Divide by LCD */
    mpz_t result;
    mpz_init(result);
    mpz_tdiv_q(result, total_sum, LCD);

    /* Verify exact division */
    mpz_t remainder;
    mpz_init(remainder);
    mpz_tdiv_r(remainder, total_sum, LCD);
    if (mpz_sgn(remainder) != 0) {
        fprintf(stderr, "ERROR: result not exactly divisible by LCD!\n");
    }
    mpz_clear(remainder);

    clock_gettime(CLOCK_MONOTONIC, &t_end);
    double elapsed = (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_nsec - t_start.tv_nsec) * 1e-9;

    /* Output */
    size_t ndigits = mpz_sizeinbase(result, 10);
    printf("a(%d) = ", n);
    mpz_out_str(stdout, 10, result);
    printf("\n");
    printf("(%zu digits, %.1fs, %lld partitions)\n", ndigits, elapsed, part_count);

    /* Cleanup */
    mpz_clear(result);
    mpz_clear(total_sum);
    mpz_clear(LCD);
    mpz_clear(tmp_z_val);
    mpz_clear(tmp_contrib);
    mpz_clear(tmp_lcd_over_z);
    mpz_clear(tmp_pow2);
    mpz_clear(tmp_factor);
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
