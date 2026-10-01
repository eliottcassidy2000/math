// opus S15 (2026-10-01), helper for shaved_tournaments_20261001.py (check T5b): which e-arc forward
// subgraphs S of TT8 (arcs i->j, i<j) embed in every 8-tournament?  Usage: u8 e classes8.txt [outfile]
// (prints the number of unavoidable S over all C(28,e) forward labellings; writes their 28-bit masks).
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#define N 8
static int nc; static uint8_t TO[7000][N];
static int pi_[28], pj_[28];
static int embed(const uint8_t *pred, const uint8_t *out){
    // map S-vertex k (in order 0..7) to T-vertex img[k]; constraint: img[j] -> img[k] for j in pred[k]
    int img[N]; uint32_t candstack[N]; int k=0; uint32_t used=0;
    uint32_t cand = 0xFF;
    candstack[0]=cand;
    while(1){
        if(candstack[k]==0){ if(k==0) return 0; k--; used &= ~(1u<<img[k]); continue; }
        int y=__builtin_ctz(candstack[k]); candstack[k]&=candstack[k]-1;
        img[k]=y; used|=1u<<y;
        if(k==N-1) return 1;
        k++;
        uint32_t c = 0xFF & ~used; uint8_t p=pred[k];
        while(p){ int j=__builtin_ctz(p); p&=p-1; c &= out[img[j]]; }
        candstack[k]=c;
    }
}
int main(int argc,char**argv){
    int e=atoi(argv[1]); FILE*f=fopen(argv[2],"r"); FILE*fo= argc>3? fopen(argv[3],"w"):NULL;
    unsigned long long code;
    while(fscanf(f,"%llu",&code)==1){ int b=0; memset(TO[nc],0,N);
        for(int i=0;i<N;i++) for(int j=i+1;j<N;j++){ if(code>>b&1) TO[nc][i]|=1<<j; else TO[nc][j]|=1<<i; b++; }
        nc++; }
    int t=0; for(int i=0;i<N;i++) for(int j=i+1;j<N;j++){ pi_[t]=i; pj_[t]=j; t++; }
    long long total=0, good=0;
    // balanced work list: (lead, second) = the two highest set bits; enumerate the remaining e-2 bits below 'second'
    int pairs[400][2], np=0;
    for(int lead=27; lead>=1; lead--) for(int sec=lead-1; sec>=0; sec--) if(sec >= e-2) { pairs[np][0]=lead; pairs[np][1]=sec; np++; }
    #pragma omp parallel reduction(+:total,good)
    {
        int ord[7000]; for(int i=0;i<nc;i++) ord[i]=i;
        #pragma omp for schedule(dynamic,1)
        for(int w=0; w<np; w++){
            int lead=pairs[w][0], sec=pairs[w][1];
            int r=e-2;
            uint32_t m = r>0 ? ((1u<<r)-1) : 0;
            while(1){
                uint32_t S = m | (1u<<lead) | (1u<<sec);
                total++;
                uint8_t pred[N]={0};
                for(uint32_t x=S; x; x&=x-1){ int b=__builtin_ctz(x); pred[pj_[b]] |= 1<<pi_[b]; }
                int ok=1;
                for(int q=0;q<nc;q++){ int ci=ord[q];
                    if(!embed(pred,TO[ci])){ ok=0; if(q){ int tmp=ord[q]; memmove(ord+1,ord,q*sizeof(int)); ord[0]=tmp; } break; } }
                if(ok){ good++;
                    #pragma omp critical
                    { if(fo) fprintf(fo,"%u\n",S); }
                }
                if(r==0) break;
                uint32_t c = m & -m, rr = m + c; m = (((rr ^ m) >> 2) / c) | rr;
                if(m >> sec) break;
            }
        }
    }
    printf("e=%d forward subsets=%lld unavoidable=%lld\n", e, total, good);
    if(fo) fclose(fo);
    return 0;
}
