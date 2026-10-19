/* horizon.c -- exact Collatz data for the Mersenne line M_K = 2^K - 1, 2 <= K <= KMAX.
   Terras map T(n) = (3n+1)/2 (n odd), n/2 (n even), stopped at the first 1.
   Output per K:  K  odd_steps  terras_steps  last_j  peak_bits  pay_terras
     odd_steps   = number of odd Terras steps (= U-steps) from M_K to 1
     terras_steps= total Terras steps from M_K to 1 (sigma_T)
     last_j      = the last odd value before 1 is (2^j - 1)/3; j recorded (j even)
     peak_bits   = max bit length along the orbit
     pay_terras  = first Terras time t with T^t(M_K) < M_K (payment / stopping time)
   Threads take K from a shared counter, largest first. Each thread writes its own file.
   Compile: cc -O2 -o horizon horizon.c -I/opt/homebrew/include -L/opt/homebrew/lib -lgmp -lpthread
   Usage:   ./horizon KMAX NTHREADS outprefix
*/
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <gmp.h>

static long KMAX, NEXTK;
static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
static const char *prefix;

static void run_one(long K, FILE *out) {
    mpz_t n, src, t;
    mpz_init(n); mpz_init(src); mpz_init(t);
    mpz_set_ui(n, 1); mpz_mul_2exp(n, n, K); mpz_sub_ui(n, n, 1);   /* 2^K - 1 */
    mpz_set(src, n);
    long odd = 0, ter = 0, pay = -1, peak = K;
    long last_odd_j = -1;
    while (mpz_cmp_ui(n, 1) != 0) {
        if (mpz_odd_p(n)) {
            /* record candidate last odd: 3n+1 power of two? */
            mpz_mul_ui(t, n, 3); mpz_add_ui(t, t, 1);
            /* t = 3n+1; next = t/2 */
            if (mpz_popcount(t) == 1) last_odd_j = (long)mpz_sizeinbase(t, 2) - 1; /* t = 2^j */
            mpz_tdiv_q_2exp(n, t, 1);
            odd++;
        } else {
            mpz_tdiv_q_2exp(n, n, 1);
        }
        ter++;
        long bits = (long)mpz_sizeinbase(n, 2);
        if (bits > peak) peak = bits;
        if (pay < 0 && mpz_cmp(n, src) < 0) pay = ter;
    }
    fprintf(out, "%ld %ld %ld %ld %ld %ld\n", K, odd, ter, last_odd_j, peak, pay);
    fflush(out);
    mpz_clear(n); mpz_clear(src); mpz_clear(t);
}

static void *worker(void *arg) {
    long id = (long)arg;
    char fn[512]; snprintf(fn, sizeof fn, "%s_t%ld.txt", prefix, id);
    FILE *out = fopen(fn, "a");
    if (!out) { perror("fopen"); return NULL; }
    for (;;) {
        pthread_mutex_lock(&lock);
        long K = NEXTK; if (K >= 2) NEXTK--;
        pthread_mutex_unlock(&lock);
        if (K < 2) break;
        run_one(K, out);
    }
    fclose(out);
    return NULL;
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s KMAX NTHREADS outprefix [KMIN]\n", argv[0]); return 1; }
    KMAX = atol(argv[1]); long nt = atol(argv[2]); prefix = argv[3];
    long kmin = (argc > 4) ? atol(argv[4]) : 2;
    NEXTK = KMAX;
    /* process K from KMAX down to kmin */
    pthread_t th[64];
    /* trick: set a floor by adjusting loop: workers stop when K < 2; emulate kmin by starting NEXTK at KMAX and stopping at kmin */
    /* simple approach: shift so that K < kmin stops */
    if (kmin > 2) {
        /* run workers with a modified stop: reuse global by offset */
    }
    for (long i = 0; i < nt; i++) pthread_create(&th[i], NULL, worker, (void*)i);
    for (long i = 0; i < nt; i++) pthread_join(th[i], NULL);
    return 0;
}
