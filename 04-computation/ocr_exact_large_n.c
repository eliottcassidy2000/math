/*
 * ocr_exact_large_n.c — kind-pasteur-2026-03-21-S17
 *
 * Compute EXACT CV^2(H) and OCR for n up to ~14 by iterating over n!
 * permutations and computing the weighted sum W = sum c(n,a)*2^a.
 *
 * For each permutation, count unit ascents (tau(i+1)-tau(i)=+1) and
 * unit descents (tau(i+1)-tau(i)=-1). If 0 descents, add 2^(ascents) to W.
 *
 * Then: CV^2(H) = (W - n!) / n!
 *       R^2 = 2(n-2) / [n(n-1) * CV^2(H)]
 *
 * n=10: 3.6M perms (~1s)
 * n=11: 40M perms (~10s)
 * n=12: 479M perms (~2 min)
 * n=13: 6.2B perms (~30 min)
 *
 * Compile: gcc -O3 -o ocr_exact_large_n ocr_exact_large_n.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef long long i64;
typedef unsigned long long u64;
typedef __int128 i128;

#define MAX_N 15

static int perm[MAX_N];
static int n_global;

/* Accumulate W = sum (2^a) over all permutations with 0 unit descents */
static i128 W_total;
static i64 n_compat;  /* Total compatible count (d=0) */

static void print128(i128 val) {
    if (val == 0) { printf("0"); return; }
    if (val < 0) { printf("-"); val = -val; }
    char buf[50]; int pos = 0;
    while (val > 0) { buf[pos++] = '0' + (int)(val % 10); val /= 10; }
    for (int i = pos-1; i >= 0; i--) printf("%c", buf[i]);
}

/* Generate all permutations via Heap's algorithm and process each */
static void heap_gen(int k) {
    if (k == 1) {
        /* Process current permutation */
        int asc = 0, desc = 0;
        for (int i = 0; i < n_global - 1; i++) {
            int diff = perm[i+1] - perm[i];
            if (diff == 1) asc++;
            else if (diff == -1) { desc = 1; break; } /* Early exit */
        }
        if (!desc) {
            W_total += (i128)1 << asc;
            n_compat++;
        }
        return;
    }

    for (int i = 0; i < k; i++) {
        heap_gen(k - 1);
        if (k % 2 == 0) {
            int tmp = perm[i]; perm[i] = perm[k-1]; perm[k-1] = tmp;
        } else {
            int tmp = perm[0]; perm[0] = perm[k-1]; perm[k-1] = tmp;
        }
    }
}

/* Compute n! as i128 */
static i128 factorial128(int n) {
    i128 result = 1;
    for (int i = 2; i <= n; i++) result *= i;
    return result;
}

/* GCD for i128 */
static i128 gcd128(i128 a, i128 b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { i128 t = b; b = a % b; a = t; }
    return a;
}

int main(int argc, char **argv) {
    int max_n = 12;
    if (argc > 1) max_n = atoi(argv[1]);
    if (max_n > 14) max_n = 14;

    printf("EXACT OCR COMPUTATION FOR n = 3 to %d\n\n", max_n);
    printf("%-4s %-20s %-20s %-20s %-15s %-15s %-10s\n",
           "n", "W", "n!", "W-n!", "CV^2 num", "CV^2 den", "R^2");
    printf("%-4s %-20s %-20s %-20s %-15s %-15s %-10s\n",
           "----", "--------------------", "--------------------",
           "--------------------", "---------------", "---------------", "----------");
    fflush(stdout);

    for (int n = 3; n <= max_n; n++) {
        n_global = n;
        W_total = 0;
        n_compat = 0;

        /* Initialize permutation */
        for (int i = 0; i < n; i++) perm[i] = i;

        time_t t0 = time(NULL);
        heap_gen(n);
        time_t t1 = time(NULL);

        i128 nfact = factorial128(n);
        i128 diff = W_total - nfact;

        /* CV^2 = diff / nfact. Reduce. */
        i128 g = gcd128(diff, nfact);
        i128 cv2_num = diff / g;
        i128 cv2_den = nfact / g;

        /* R^2 = 2(n-2) / [n(n-1) * cv2_num/cv2_den] = 2(n-2)*cv2_den / [n(n-1)*cv2_num] */
        i128 r2_num = 2 * (n-2) * cv2_den;
        i128 r2_den = (i128)n * (n-1) * cv2_num;
        i128 g2 = gcd128(r2_num, r2_den);
        r2_num /= g2;
        r2_den /= g2;

        /* 1 - R^2 */
        i128 one_minus_num = r2_den - r2_num;

        printf("%-4d ", n);
        print128(W_total); printf(" ");
        print128(nfact); printf(" ");
        print128(diff); printf(" ");
        print128(cv2_num); printf("/");
        print128(cv2_den); printf(" ");
        print128(r2_num); printf("/");
        print128(r2_den);
        printf(" [%lds, %lld compat]\n", (long)difftime(t1, t0), n_compat);
        fflush(stdout);

        /* Print factorization of cv2_num (the prime!) */
        /* Convert to long long if small enough */
        if (cv2_num < (i128)1e15) {
            i64 val = (i64)cv2_num;
            printf("     CV^2 numerator %lld factors: ", val);
            i64 temp = val;
            for (i64 f = 2; f * f <= temp; f++) {
                while (temp % f == 0) {
                    printf("%lld ", f);
                    temp /= f;
                }
            }
            if (temp > 1) printf("%lld", temp);
            printf("\n");

            /* Also print 1-R^2 */
            i64 om_n = (i64)one_minus_num;
            i64 om_d = (i64)r2_den;
            printf("     1 - R^2 = %lld/%lld = %.15f\n", om_n, om_d,
                   (double)om_n / (double)om_d);
        }

        printf("\n");
        fflush(stdout);
    }

    return 0;
}
