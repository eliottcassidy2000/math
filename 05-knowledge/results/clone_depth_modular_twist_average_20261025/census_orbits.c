/* Exact census of orbits of a regular order-25 group action on subsets of its 25 points, by subset size.
 * Independent of Burnside: an orbit is counted once, at its lexicographically minimal mask.
 * Usage: ./census_orbits z5sq   (translations of Z5 x Z5)   |   ./census_orbits z25   (rotations of Z25)
 * Build: cc -O2 -o census_orbits census_orbits.c
 */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#define P 5
#define N 25
static int perm[N][N];
static uint32_t tab[N][4][256];
int main(int argc, char **argv){
    int cyclic = (argc > 1 && strcmp(argv[1], "z25") == 0);
    for(int g=0; g<N; g++) for(int i=0;i<N;i++){
        if(cyclic) perm[g][i] = (i + g) % N;
        else { int gc=g/P, gd=g%P, a=i/P, b=i%P; perm[g][i]=((a+gc)%P)*P+((b+gd)%P); }
    }
    for(int g=0; g<N; g++) for(int byte=0; byte<4; byte++) for(int v=0; v<256; v++){
        uint32_t m=0; for(int bit=0;bit<8;bit++) if(v>>bit&1){ int i=byte*8+bit; if(i<N) m|=1u<<perm[g][i]; } tab[g][byte][v]=m; }
    unsigned long long count[N+1]; memset(count,0,sizeof count); unsigned long long orbits=0;
    for(uint32_t mask=0; mask < (1u<<N); mask++){
        int minimal=1;
        for(int g=1; g<N; g++){
            uint32_t im = tab[g][0][mask&255] | tab[g][1][(mask>>8)&255] | tab[g][2][(mask>>16)&255] | tab[g][3][(mask>>24)&255];
            if(im < mask){ minimal=0; break; }
        }
        if(minimal){ count[__builtin_popcount(mask)]++; orbits++; }
    }
    printf("action %s\norbits %llu\ncoeffs", cyclic ? "Z25 rotations" : "Z5xZ5 translations", orbits);
    for(int k=0;k<=N;k++) printf(" %llu", count[k]);
    printf("\n");
    return 0;
}
