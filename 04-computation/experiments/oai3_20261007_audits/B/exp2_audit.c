// Independent audit: negative discriminants D = -N (N = 0,3 mod 4) whose form class group
// (primitive forms) has exponent <= 2.  Exponent <= 2 <=> no primitive reduced form (a,b,c)
// with 0 < b < a < c (such a class differs from its inverse (a,-b,c)).
// Part 1: direct early-exit check for N <= N1.
// Part 2: sieve for N1 < N <= N2: if p odd prime, p !| N, (-N/p) = 1 and 4p^2 < N, the form
// (p, b, (b^2+N)/(4p)) with 0<|b|<p is primitive reduced non-ambiguous; also N = 7 mod 8 (D = 1 mod 8)
// with N > 15 gives (2,1,c). Survivors get the direct check.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

static uint64_t gcd(uint64_t a, uint64_t b){ while(b){ uint64_t t=a%b; a=b; b=t;} return a; }

// returns 1 if exponent <= 2 (no non-ambiguous primitive reduced form)
static int exp_le2(uint64_t N){
    for (uint64_t a = 1; 3*a*a <= N; a++){
        uint64_t four_a = 4*a;
        for (uint64_t b = (N & 1) ? 1 : 2; b < a; b += 2){
            uint64_t t = b*b + N;
            if (t % four_a) continue;
            uint64_t c = t / four_a;
            if (c <= a) continue;
            if (gcd(gcd(a,b),c) != 1) continue;
            return 0;
        }
    }
    return 1;
}

int main(int argc, char **argv){
    uint64_t N1 = strtoull(argv[1],0,10), N2 = strtoull(argv[2],0,10);
    int found = 0;
    for (uint64_t N = 3; N <= N1; N++){
        if ((N & 3) != 0 && (N & 3) != 3) continue;
        if (exp_le2(N)){ printf("%llu\n",(unsigned long long)N); found++; }
    }
    fprintf(stderr,"direct part N<=%llu: %d\n",(unsigned long long)N1,found);
    fflush(stdout);
    if (N2 <= N1) return 0;
    // sieve part
    int P = 400; int primes[200]; int np=0;
    for (int p=3;p<=P;p++){ int ok=1; for(int d=2;d*d<=p;d++) if(p%d==0){ok=0;break;} if(ok) primes[np++]=p; }
    // bad[p][r] = 1 if r = N mod p makes (-N/p) = 1
    static unsigned char bad[400][400];
    for (int i=0;i<np;i++){ int p=primes[i]; memset(bad[p],0,400);
        for (int x=1;x<p;x++){ int s=(x*x)%p; int r=(p - s)%p; bad[p][r]=1; } }  // -N = s (QR) <=> N = -s
    // wheel modulus over 8 and primes 3..17
    const uint64_t M = 8ULL*3*5*7*11*13*17;
    uint32_t *res = malloc(sizeof(uint32_t)*M/8); int nres=0;
    for (uint64_t r=0;r<M;r++){
        int m8 = r & 7; if (!(m8==0||m8==3||m8==4)) continue;
        int ok=1; for (int i=0;i<6;i++){ int p=primes[i]; if (bad[p][r%p]) {ok=0;break;} }
        if (ok) res[nres++]=(uint32_t)r;
    }
    fprintf(stderr,"wheel residues %d of %llu\n",nres,(unsigned long long)M);
    uint64_t cand=0, surv=0; int found2=0;
    uint64_t q0 = (N1+1)/M, q1 = N2/M;
    for (uint64_t q=q0; q<=q1; q++){
        uint64_t base=q*M;
        // residues of base mod each larger prime
        int baser[200]; for (int i=6;i<np;i++) baser[i]=(int)(base % primes[i]);
        for (int k=0;k<nres;k++){
            uint64_t N = base + res[k];
            if (N <= N1 || N > N2) continue;
            cand++;
            int ok=1;
            for (int i=6;i<np;i++){ int p=primes[i]; if ((uint64_t)4*p*p >= N) break;
                int r = (baser[i] + (int)(res[k] % p)) % p; if (bad[p][r]){ok=0;break;} }
            if (!ok) continue;
            surv++;
            if (exp_le2(N)){ printf("%llu\n",(unsigned long long)N); found2++; fflush(stdout);}
        }
    }
    fprintf(stderr,"sieve part (%llu,%llu]: candidates %llu survivors %llu found %d\n",(unsigned long long)N1,(unsigned long long)N2,(unsigned long long)cand,(unsigned long long)surv,found2);
    return 0;
}
