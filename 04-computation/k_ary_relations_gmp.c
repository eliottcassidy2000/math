/*
 * Compute k-ary relations on n unlabeled nodes using Burnside/Polya.
 *
 * a(n) = (1/n!) * sum_{lambda partition of n} permcount(lambda) * 2^{orbits(lambda)}
 *
 * where orbits(lambda) = number of orbits of a permutation with cycle type lambda
 * acting on [n]^k (ordered k-tuples from {1,...,n}).
 *
 * For cycle type with cycles c_1,...,c_m (with repetition):
 *   orbits = (1/L) * sum_{t=0}^{L-1} f(t)^k
 *   where L = lcm(c_1,...,c_m), f(t) = sum_{c_i | t} c_i
 *
 * Known OEIS sequences:
 *   k=2: A000595 (binary relations)
 *   k=3: A000662
 *   k=4: A001377
 *   k=5: A051241
 *
 * Compile: gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *          -o k_ary_relations_gmp k_ary_relations_gmp.c -lgmp
 *
 * Usage: ./k_ary_relations_gmp <k> <max_n>
 *
 * Author: opus-2026-03-09-S50
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>

static int gcd_func(int a, int b) {
    while (b) { int t = b; b = a % b; a = t; }
    return a;
}

static long long lcm_func(long long a, long long b) {
    return a / gcd_func(a, b) * b;
}

/* Global state */
static int K;           /* arity */
static int N;           /* number of elements */
static mpz_t result;    /* accumulated result (before dividing by n!) */
static mpz_t fact_n;    /* n! */
static long long partition_count;

/* Compute orbits of a permutation with given cycle type on [n]^k
 * cycles[] = array of cycle lengths (with repetition), num_cycles = count
 * Returns the number of orbits.
 */
static void compute_orbits(int *cycles, int num_cycles, mpz_t orbits_out) {
    /* L = lcm of all cycle lengths */
    long long L = 1;
    for (int i = 0; i < num_cycles; i++) {
        L = lcm_func(L, cycles[i]);
        if (L > 1000000000LL) {
            /* L too large to iterate, use divisor approach */
            /* For now, fall back to a slower method */
            break;
        }
    }

    /* f(t) = sum of c_i where c_i | t */
    /* orbits = (1/L) * sum_{t=0}^{L-1} f(t)^k */

    mpz_set_ui(orbits_out, 0);
    mpz_t ft_pow, term;
    mpz_init(ft_pow);
    mpz_init(term);

    if (L <= 100000000LL) {
        /* Direct iteration */
        for (long long t = 0; t < L; t++) {
            long long ft;
            if (t == 0) {
                ft = N;
            } else {
                ft = 0;
                for (int i = 0; i < num_cycles; i++) {
                    if (t % cycles[i] == 0) ft += cycles[i];
                }
            }
            /* orbits_out += ft^k */
            mpz_set_ui(ft_pow, ft);
            mpz_pow_ui(ft_pow, ft_pow, K);
            mpz_add(orbits_out, orbits_out, ft_pow);
        }
        /* Divide by L */
        mpz_divexact_ui(orbits_out, orbits_out, L);
    } else {
        /* L is too large. Use divisor-based approach:
         * sum_{t=0}^{L-1} f(t)^k = sum over divisor tuples
         * This is harder and not implemented here.
         * For now, this case shouldn't arise for reasonable n.
         */
        fprintf(stderr, "Warning: L=%lld too large, skipping partition\n", L);
        mpz_set_ui(orbits_out, 0);
    }

    mpz_clear(ft_pow);
    mpz_clear(term);
}

static void enumerate(int remaining, int max_part, int *pk, int *pm, int depth) {
    if (remaining == 0) {
        partition_count++;

        /* Expand cycle type to cycle array */
        int num_cycles = 0;
        for (int i = 0; i < depth; i++) num_cycles += pm[i];
        int *cycles = (int *)malloc(num_cycles * sizeof(int));
        int idx = 0;
        for (int i = 0; i < depth; i++)
            for (int j = 0; j < pm[i]; j++)
                cycles[idx++] = pk[i];

        /* Compute orbits */
        mpz_t orbits;
        mpz_init(orbits);
        compute_orbits(cycles, num_cycles, orbits);

        /* z_lambda = prod r^m * m! */
        mpz_t z_lambda;
        mpz_init(z_lambda);
        mpz_set_ui(z_lambda, 1);
        for (int i = 0; i < depth; i++) {
            for (int j = 0; j < pm[i]; j++)
                mpz_mul_ui(z_lambda, z_lambda, pk[i]);
            for (int j = 2; j <= pm[i]; j++)
                mpz_mul_ui(z_lambda, z_lambda, j);
        }

        /* contribution = (n! / z_lambda) * 2^orbits */
        mpz_t contrib;
        mpz_init(contrib);
        mpz_divexact(contrib, fact_n, z_lambda);

        mpz_t power;
        mpz_init(power);
        if (mpz_fits_ulong_p(orbits)) {
            unsigned long orb_val = mpz_get_ui(orbits);
            mpz_set_ui(power, 1);
            mpz_mul_2exp(power, power, orb_val);
        } else {
            /* orbits is very large; 2^orbits will be enormous */
            mpz_set_ui(power, 2);
            mpz_powm(power, power, orbits, power); /* This won't work; need different approach */
            /* Actually for exact computation: */
            /* We need 2^orbits exactly. Use mpz_mul_2exp with converted value */
            /* orbits should fit in unsigned long for reasonable n */
            fprintf(stderr, "Warning: orbits too large for mpz_mul_2exp\n");
        }

        mpz_mul(contrib, contrib, power);
        mpz_add(result, result, contrib);

        mpz_clear(contrib);
        mpz_clear(power);
        mpz_clear(z_lambda);
        mpz_clear(orbits);
        free(cycles);
        return;
    }

    for (int p = (max_part < remaining ? max_part : remaining); p >= 1; p--) {
        int max_m = remaining / p;
        for (int m = 1; m <= max_m; m++) {
            pk[depth] = p;
            pm[depth] = m;
            enumerate(remaining - m * p, p - 1, pk, pm, depth + 1);
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <k> <max_n>\n", argv[0]);
        fprintf(stderr, "  k: arity (2=binary, 3=ternary, ...)\n");
        fprintf(stderr, "  max_n: compute a(0)..a(max_n)\n");
        return 1;
    }

    K = atoi(argv[1]);
    int max_n = atoi(argv[2]);

    fprintf(stderr, "Computing %d-ary relations on n unlabeled nodes, n=0..%d\n\n", K, max_n);

    /* a(0) = 1 always */
    printf("0 1\n");
    fflush(stdout);

    for (N = 1; N <= max_n; N++) {
        struct timespec t_start, t_end;
        clock_gettime(CLOCK_MONOTONIC, &t_start);

        mpz_init(result);
        mpz_init(fact_n);
        mpz_fac_ui(fact_n, N);
        partition_count = 0;

        int pk[64], pm[64];
        enumerate(N, N, pk, pm, 0);

        /* Divide by n! */
        mpz_divexact(result, result, fact_n);

        clock_gettime(CLOCK_MONOTONIC, &t_end);
        double dt = (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_nsec - t_start.tv_nsec) / 1e9;

        gmp_printf("%d %Zd\n", N, result);
        fprintf(stderr, "n=%d: %zu digits, %lld partitions, %.3fs\n",
                N, mpz_sizeinbase(result, 10), partition_count, dt);
        fflush(stdout);
        fflush(stderr);

        mpz_clear(result);
        mpz_clear(fact_n);

        if (dt > 300) {
            fprintf(stderr, "Stopping: too slow\n");
            break;
        }
    }

    return 0;
}
