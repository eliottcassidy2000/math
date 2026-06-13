/*
 * A008406: Triangle T(n,k) = number of simple graphs on n nodes with k edges.
 *
 * Uses Burnside/Polya: T(n,k) = [x^k] (1/n!) sum_{sigma in S_n} GF_sigma(x)
 *
 * For cycle type lambda = [(r_a, m_a)]:
 *   GF_sigma(x) = prod of (1 + x^{orbit_size}) over all pair orbits
 *
 * Orbit structure:
 *   Within single r-cycle: (r-1)/2 orbits of size r (r odd)
 *                          (r-2)/2 orbits of size r + 1 orbit of size r/2 (r even)
 *   Within m copies of r-cycle: C(m,2)*r orbits of size r
 *   Cross (r_a, m_a) and (r_b, m_b): m_a*m_b*gcd(r_a,r_b) orbits of size lcm(r_a,r_b)
 *
 * Polynomial arithmetic in GMP for exact results.
 * Output: b-file format for OEIS.
 *
 * Compile: gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *          -o a008406_gmp a008406_gmp.c -lgmp
 *
 * Usage: ./a008406_gmp <max_n>
 *
 * Author: opus-2026-03-08-S48
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>

static int gcd(int a, int b) {
    while (b) { int t = b; b = a % b; a = t; }
    return a;
}

/*
 * Multiply polynomial poly[0..max_deg] by (1 + x^d)^e in-place.
 * Uses the identity: multiplying by (1+x^d) shifts and adds.
 * For efficiency, we iterate e times, each time doing poly *= (1+x^d).
 */
static void poly_mul_1_plus_xd(mpz_t *poly, int max_deg, int d, int e) {
    if (e == 0 || d > max_deg) return;

    /* For small e, iterate: poly *= (1+x^d) e times */
    /* For large e, expand (1+x^d)^e as binomial and convolve */

    if (e <= 64 || (long long)e * d <= 4 * max_deg) {
        /* Iterative: for each factor (1+x^d), update poly from high to low */
        for (int iter = 0; iter < e; iter++) {
            for (int i = max_deg; i >= d; i--) {
                mpz_add(poly[i], poly[i], poly[i - d]);
            }
        }
    } else {
        /* For very large e: use binomial expansion */
        /* (1+x^d)^e = sum_{j=0}^{e} C(e,j) x^{j*d} */
        int max_j = max_deg / d;
        if (max_j > e) max_j = e;

        /* Build binomial coefficients */
        mpz_t *binom = (mpz_t *)malloc((max_j + 1) * sizeof(mpz_t));
        for (int j = 0; j <= max_j; j++) mpz_init(binom[j]);
        mpz_set_ui(binom[0], 1);
        for (int j = 1; j <= max_j; j++) {
            mpz_set(binom[j], binom[j-1]);
            mpz_mul_ui(binom[j], binom[j], e - j + 1);
            mpz_divexact_ui(binom[j], binom[j], j);
        }

        /* Convolve: result[i] = sum_{j} poly_old[i - j*d] * binom[j] */
        mpz_t *result = (mpz_t *)malloc((max_deg + 1) * sizeof(mpz_t));
        for (int i = 0; i <= max_deg; i++) mpz_init_set_ui(result[i], 0);

        mpz_t tmp;
        mpz_init(tmp);
        for (int i = 0; i <= max_deg; i++) {
            if (mpz_sgn(poly[i]) == 0) continue;
            for (int j = 0; j <= max_j; j++) {
                int idx = i + j * d;
                if (idx > max_deg) break;
                mpz_mul(tmp, poly[i], binom[j]);
                mpz_add(result[idx], result[idx], tmp);
            }
        }
        mpz_clear(tmp);

        /* Copy result back */
        for (int i = 0; i <= max_deg; i++) {
            mpz_set(poly[i], result[i]);
            mpz_clear(result[i]);
        }
        free(result);

        for (int j = 0; j <= max_j; j++) mpz_clear(binom[j]);
        free(binom);
    }
}

/*
 * Compute GF(x) for a partition in compressed form.
 * compressed[i] = (r, m), depth = number of distinct parts.
 */
static void compute_gf(int *pk, int *pm, int depth, mpz_t *poly, int max_deg) {
    /* Initialize polynomial to 1 */
    for (int i = 0; i <= max_deg; i++) mpz_set_ui(poly[i], 0);
    mpz_set_ui(poly[0], 1);

    for (int idx = 0; idx < depth; idx++) {
        int r = pk[idx], m = pm[idx];

        /* Within-cycle terms */
        if (r % 2 == 1) {
            /* r odd: (1+x^r)^{m*(r-1)/2} */
            int exp = m * (r - 1) / 2;
            if (exp > 0)
                poly_mul_1_plus_xd(poly, max_deg, r, exp);
        } else {
            /* r even: (1+x^r)^{m*(r-2)/2} * (1+x^{r/2})^m */
            int exp1 = m * (r - 2) / 2;
            if (exp1 > 0)
                poly_mul_1_plus_xd(poly, max_deg, r, exp1);
            poly_mul_1_plus_xd(poly, max_deg, r / 2, m);
        }

        /* Same-part cross-cycle terms: C(m,2)*r orbits of size r */
        int cross_same = m * (m - 1) / 2 * r;
        if (cross_same > 0)
            poly_mul_1_plus_xd(poly, max_deg, r, cross_same);
    }

    /* Cross-part terms */
    for (int i = 0; i < depth; i++) {
        int ri = pk[i], mi = pm[i];
        for (int j = i + 1; j < depth; j++) {
            int rj = pk[j], mj = pm[j];
            int g = gcd(ri, rj);
            int l = ri / g * rj;  /* lcm */
            int exp = mi * mj * g;
            if (exp > 0 && l <= max_deg)
                poly_mul_1_plus_xd(poly, max_deg, l, exp);
        }
    }
}

/* Global state for partition enumeration */
static int global_n, global_max_edges;
static mpz_t *total_poly;  /* accumulated weighted GF */
static mpz_t z_lambda;     /* automorphism group order */
static long long partition_count;

static void enumerate(int remaining, int max_part, int *pk, int *pm, int depth) {
    if (remaining == 0) {
        partition_count++;

        /* Compute GF for this partition */
        mpz_t *gf = (mpz_t *)malloc((global_max_edges + 1) * sizeof(mpz_t));
        for (int i = 0; i <= global_max_edges; i++) mpz_init(gf[i]);

        compute_gf(pk, pm, depth, gf, global_max_edges);

        /* Compute z_lambda = prod r^m * m! */
        mpz_set_ui(z_lambda, 1);
        for (int i = 0; i < depth; i++) {
            for (int j = 0; j < pm[i]; j++)
                mpz_mul_ui(z_lambda, z_lambda, pk[i]);
            for (int j = 2; j <= pm[i]; j++)
                mpz_mul_ui(z_lambda, z_lambda, j);
        }

        /* Accumulate: total_poly[k] += gf[k] / z_lambda */
        /* We work with common denominator later; for now accumulate as rationals */
        /* Actually, let's use mpq for exact arithmetic */
        /* Better approach: accumulate numerators with LCD */

        /* For simplicity, add gf[k] * (n!/z_lambda) to total, then divide by n! at end */
        mpz_t scale;
        mpz_init(scale);
        mpz_fac_ui(scale, global_n);
        mpz_divexact(scale, scale, z_lambda);

        for (int k = 0; k <= global_max_edges; k++) {
            if (mpz_sgn(gf[k]) != 0) {
                mpz_addmul(total_poly[k], gf[k], scale);
            }
        }

        mpz_clear(scale);
        for (int i = 0; i <= global_max_edges; i++) mpz_clear(gf[i]);
        free(gf);
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
    int max_n = argc > 1 ? atoi(argv[1]) : 20;

    printf("A008406: Triangle T(n,k) = graphs on n nodes with k edges\n");
    printf("Computing rows 1..%d\n\n", max_n);

    /* Open b-file for output */
    FILE *bfile = fopen("b008406.txt", "w");
    int bfile_idx = 0;

    for (int n = 1; n <= max_n; n++) {
        global_n = n;
        global_max_edges = n * (n - 1) / 2;
        partition_count = 0;

        /* Allocate total polynomial */
        total_poly = (mpz_t *)malloc((global_max_edges + 1) * sizeof(mpz_t));
        for (int i = 0; i <= global_max_edges; i++)
            mpz_init_set_ui(total_poly[i], 0);
        mpz_init(z_lambda);

        int pk[64], pm[64];

        clock_t t0 = clock();
        enumerate(n, n, pk, pm, 0);
        double dt = (double)(clock() - t0) / CLOCKS_PER_SEC;

        /* Divide by n! */
        mpz_t fact_n;
        mpz_init(fact_n);
        mpz_fac_ui(fact_n, n);

        /* Divide each coefficient by n! and write to b-file */
        for (int k = 0; k <= global_max_edges; k++) {
            mpz_divexact(total_poly[k], total_poly[k], fact_n);
            if (bfile) {
                fprintf(bfile, "%d ", bfile_idx);
                mpz_out_str(bfile, 10, total_poly[k]);
                fprintf(bfile, "\n");
                bfile_idx++;
            }
        }

        /* Print summary */
        printf("n=%2d: %d terms, T(n,0)=", n, global_max_edges + 1);
        mpz_out_str(stdout, 10, total_poly[0]);
        if (global_max_edges > 0) {
            printf(", T(n,1)=");
            mpz_out_str(stdout, 10, total_poly[1]);
        }
        if (global_max_edges > 4) {
            printf(", ..., T(n,C(n,2))=");
            mpz_out_str(stdout, 10, total_poly[global_max_edges]);
        }
        printf("  [%.3fs, %lld parts]\n", dt, partition_count);
        fflush(stdout);

        mpz_clear(fact_n);
        mpz_clear(z_lambda);
        for (int i = 0; i <= global_max_edges; i++)
            mpz_clear(total_poly[i]);
        free(total_poly);

        if (dt > 300) {
            printf("Stopping: too slow\n");
            break;
        }
    }

    if (bfile) {
        fclose(bfile);
        printf("\nWritten b008406.txt with %d entries\n", bfile_idx);
    }

    return 0;
}
