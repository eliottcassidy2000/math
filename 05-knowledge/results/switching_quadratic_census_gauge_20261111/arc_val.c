/* 2-adic valuations of mixed arc-flip derivatives of H, by shape of the flipped set F (a linear forest).
   usage: arc_val n      (n <= 7; all labelled tournaments)
   D_F H(T) = sum_{G subset F} (-1)^{|F|-|G|} H(T with the arcs of G reversed); min over ALL T of v_2(D_F H(T)). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int n, npairs, pidx[8][8];
static int32_t Hc(uint32_t T){
    int out[8]; int full=(1<<n)-1; for(int i=0;i<n;i++) out[i]=0;
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if((T>>pidx[i][j])&1) out[i]|=1<<j; else out[j]|=1<<i; }
    static int32_t dp[128][7]; memset(dp,0,sizeof(dp)); for(int v=0;v<n;v++) dp[1<<v][v]=1;
    for(int S=1;S<=full;S++) for(int v=0;v<n;v++){ int32_t c=dp[S][v]; if(!c) continue; int t=out[v]&~S; while(t){int b=t&-t; dp[S|b][__builtin_ctz(b)]+=c; t^=b;} }
    int32_t h=0; for(int v=0;v<n;v++) h+=dp[full][v]; return h;
}
static int v2(int32_t x){ if(x==0) return 99; int v=0; while(!(x&1)){x/=2;v++;} return v; }
int main(int argc,char**argv){
    n=atoi(argv[1]); npairs=0; for(int i=0;i<n;i++) for(int j=i+1;j<n;j++) pidx[i][j]=pidx[j][i]=npairs++;
    uint32_t NT=1u<<npairs; int32_t *H=malloc(NT*sizeof(int32_t)), *D=malloc(NT*sizeof(int32_t));
    for(uint32_t T=0;T<NT;T++) H[T]=Hc(T);
    /* shapes: each is a list of vertex paths */
    const char *shapes[]={"01","012","0123","01234","012345","0123456","01|23","01|23|45","012|34","012|345","0123|45","01|234|56","012|34|56","0123|456","01234|56",NULL};
    for(int s=0;shapes[s];s++){
        const char *p=shapes[s]; int bits[32]; int nb=0; int ok=1; int comp_edges[8]; int nc=0; int cur=0;
        for(int i=0;p[i];i++){
            if(p[i]=='|'){ comp_edges[nc++]=cur; cur=0; continue; }
            int a=p[i]-'0'; if(a>=n){ok=0;break;}
            if(p[i+1]&&p[i+1]!='|'){ int b=p[i+1]-'0'; if(b>=n){ok=0;break;} bits[nb++]=pidx[a][b]; cur++; }
        }
        comp_edges[nc++]=cur; if(!ok) continue;
        memcpy(D,H,NT*sizeof(int32_t));
        for(int k=0;k<nb;k++){ uint32_t b=1u<<bits[k]; for(uint32_t T=0;T<NT;T++) if(!(T&b)){ int32_t a0=D[T], a1=D[T|b]; D[T]=a1-a0; D[T|b]=a0-a1; } }
        int mn=99; long cnt=0, nz=0; for(uint32_t T=0;T<NT;T++){ int v=v2(D[T]); if(v<99) nz++; if(v<mn){mn=v;cnt=1;} else if(v==mn) cnt++; }
        int pred=0; for(int c=0;c<nc;c++) pred+=(comp_edges[c]+1)/2;
        printf("n=%d shape %-12s edges per component:",n,p); for(int c=0;c<nc;c++) printf(" %d",comp_edges[c]);
        printf("   min v_2 over all T = %d (nonzero at %ld tournaments)   sum of ceil(m_c/2) = %d   #components = %d\n",mn,nz,pred,nc);
    }
    return 0;
}
