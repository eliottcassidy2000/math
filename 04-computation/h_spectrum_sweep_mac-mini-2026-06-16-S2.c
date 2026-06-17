/* h_spectrum_sweep — mac-mini-2026-06-16-S2
 * Exhaustive realized H-spectrum (H = #directed Hamiltonian paths) over ALL
 * 2^C(n,2) labeled tournaments on n vertices, via Held-Karp DP counting.
 * Tests HYP-2558: is residue 1 mod 8 (the octal odd-square class) gap-free,
 * and do the n<=6 gaps {35,39} (and 7,21) get realized at larger n?
 * Usage: ./h_spectrum n [num_threads]
 * gcc -O3 -pthread.  H <= n! so a small bool table suffices.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#define MAXH 400000          /* > 9! = 362880 */
static int N, E;
static int pr_i[64], pr_j[64];   /* edge endpoints in combination order */
static char realized_global[MAXH];
static long maxH_global = 0;
static pthread_mutex_t mtx = PTHREAD_MUTEX_INITIALIZER;

typedef struct { long start, stop; } Range;

static void* worker(void* arg){
    Range* r = (Range*)arg;
    char* local = calloc(MAXH, 1);
    long localmax = 0;
    int full = (1<<N) - 1;
    long dpsize = (long)(1<<N) * N;
    long* dp = malloc(sizeof(long) * dpsize);
    int adj[16];
    for(long bits = r->start; bits < r->stop; bits++){
        for(int v=0; v<N; v++) adj[v]=0;
        for(int k=0; k<E; k++){
            if(bits & (1L<<k)) adj[pr_i[k]] |= (1<<pr_j[k]);
            else               adj[pr_j[k]] |= (1<<pr_i[k]);
        }
        memset(dp, 0, sizeof(long)*dpsize);
        for(int v=0; v<N; v++) dp[(long)(1<<v)*N + v] = 1;
        for(int mask=1; mask<=full; mask++){
            long base = (long)mask*N;
            for(int v=0; v<N; v++){
                if(!(mask & (1<<v))) continue;
                long c = dp[base+v];
                if(!c) continue;
                int outs = adj[v] & ~mask;
                while(outs){
                    int u = __builtin_ctz(outs);
                    outs &= outs-1;
                    dp[(long)(mask|(1<<u))*N + u] += c;
                }
            }
        }
        long H = 0;
        for(int v=0; v<N; v++) H += dp[(long)full*N + v];
        if(H < MAXH) local[H] = 1;
        if(H > localmax) localmax = H;
    }
    pthread_mutex_lock(&mtx);
    for(long h=0; h<MAXH; h++) if(local[h]) realized_global[h]=1;
    if(localmax > maxH_global) maxH_global = localmax;
    pthread_mutex_unlock(&mtx);
    free(local); free(dp);
    return NULL;
}

int main(int argc, char** argv){
    N = (argc>1)? atoi(argv[1]) : 7;
    int nt = (argc>2)? atoi(argv[2]) : 1;
    E = N*(N-1)/2;
    int k=0;
    for(int i=0;i<N;i++) for(int j=i+1;j<N;j++){ pr_i[k]=i; pr_j[k]=j; k++; }
    long total = 1L << E;
    fprintf(stderr,"n=%d E=%d tournaments=2^%d=%ld threads=%d\n", N,E,E,total,nt);

    pthread_t th[64];
    Range rg[64];
    long chunk = total / nt;
    for(int t=0;t<nt;t++){
        rg[t].start = (long)t*chunk;
        rg[t].stop  = (t==nt-1)? total : (long)(t+1)*chunk;
        pthread_create(&th[t], NULL, worker, &rg[t]);
    }
    for(int t=0;t<nt;t++) pthread_join(th[t], NULL);

    /* report */
    long maxH = maxH_global;
    printf("n=%d : max H = %ld\n", N, maxH);
    printf("realized H-values:");
    int cnt=0;
    for(long h=0;h<=maxH;h++) if(realized_global[h]){ printf(" %ld", h); cnt++; }
    printf("\n#realized = %d\n", cnt);

    printf("ODD integers in [1,%ld] NOT realized (gaps):", maxH);
    long gaps[100000]; int ng=0;
    for(long h=1;h<=maxH;h+=2) if(!realized_global[h]){ printf(" %ld", h); if(ng<100000) gaps[ng++]=h; }
    printf("\n#gaps = %d\n", ng);

    long rmod[8]={0}, gmod[8]={0};
    for(long h=1;h<=maxH;h+=2){ if(realized_global[h]) rmod[h%8]++; else gmod[h%8]++; }
    printf("realized odd by residue mod 8: 1->%ld 3->%ld 5->%ld 7->%ld\n", rmod[1],rmod[3],rmod[5],rmod[7]);
    printf("GAPS     odd by residue mod 8: 1->%ld 3->%ld 5->%ld 7->%ld\n", gmod[1],gmod[3],gmod[5],gmod[7]);
    printf("HYP-2558 (residue 1 gap-free): %s\n", gmod[1]==0 ? "HOLDS" : "BROKEN");
    /* check specific n<=6 gaps */
    long checks[4]={7,21,35,39};
    for(int c=0;c<4;c++){ long v=checks[c]; if(v<=maxH) printf("  %ld realized at n=%d? %s (res %ld)\n", v, N, realized_global[v]?"YES":"no", v%8); }
    return 0;
}
