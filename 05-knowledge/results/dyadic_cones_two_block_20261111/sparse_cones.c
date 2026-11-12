/* Cone roots (B, O): response M a power of two, any offset c0 = H(B) - W(O); derived from the unit-offset search.
   Unit-offset roots (B, O) with B = transitive tournament on 0..m-1 plus k back arcs, |O| <= 3.
   usage: sparse_roots m k [cmax]
   prints, for each odd c <= cmax, the least j with a root of response M = c 2^j, and an example. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
static int m; static int out[16], inn[16];
static u64 *E, *St, *Hs;
static u64 Hval, U[16], W[16][16];
static void compute(void){
    int N=1<<m, full=N-1;
    memset(E,0,sizeof(u64)*N*m); memset(St,0,sizeof(u64)*N*m);
    for(int v=0;v<m;v++){ E[(1<<v)*m+v]=1; St[(1<<v)*m+v]=1; }
    for(int S=1;S<N;S++){ u64*r=E+(size_t)S*m; for(int v=0;v<m;v++){ u64 c=r[v]; if(!c) continue; int t=out[v]&~S; while(t){int b=t&-t; E[(size_t)(S|b)*m+__builtin_ctz(b)]+=c; t^=b;} }
                          u64*q=St+(size_t)S*m; for(int v=0;v<m;v++){ u64 c=q[v]; if(!c) continue; int t=inn[v]&~S; while(t){int b=t&-t; St[(size_t)(S|b)*m+__builtin_ctz(b)]+=c; t^=b;} } }
    for(int S=0;S<N;S++){ u64 h=0; for(int v=0;v<m;v++) h+=E[(size_t)S*m+v]; Hs[S]=h; } Hs[0]=1;
    Hval=Hs[full]; memset(U,0,sizeof(U)); memset(W,0,sizeof(W));
    for(int S=0;S<N;S++){ int C=full^S; if(!C) continue; u64 hs=Hs[S];
        for(int y=0;y<m;y++) if((C>>y)&1) U[y]+=hs*St[(size_t)C*m+y];
        if(S) for(int e=0;e<m;e++) if((S>>e)&1){ u64 te=E[(size_t)S*m+e]; if(!te) continue; for(int f=0;f<m;f++) if((C>>f)&1) W[e][f]+=te*St[(size_t)C*m+f]; } }
}
#define CMAX 100000
static int bestj[CMAX]; static char bestdesc[CMAX][96];
int main(int argc,char**argv){
    m=atoi(argv[1]); int k=atoi(argv[2]); int cmax=argc>3?atoi(argv[3]):600; int N=1<<m;
    E=malloc(sizeof(u64)*N*m); St=malloc(sizeof(u64)*N*m); Hs=malloc(sizeof(u64)*N);
    for(int c=0;c<CMAX;c++) bestj[c]=99;
    int np=0, pa[128], pb[128]; for(int i=0;i<m;i++) for(int j=i+1;j<m;j++){ pa[np]=j; pb[np]=i; np++; }
    int idx[8]; for(int i=0;i<k;i++) idx[i]=i; long nB=0, nroots=0;
    while(1){
        for(int i=0;i<m;i++){ out[i]=0; inn[i]=0; }
        for(int i=0;i<m;i++) for(int j=i+1;j<m;j++){ out[i]|=1<<j; inn[j]|=1<<i; }
        for(int t=0;t<k;t++){ int j=pa[idx[t]], i=pb[idx[t]]; out[i]&=~(1<<j); inn[j]&=~(1<<i); out[j]|=1<<i; inn[i]|=1<<j; }
        compute(); nB++;
        {   char buf[64]; int p=0; buf[0]=0; for(int t=0;t<k;t++) p+=snprintf(buf+p,64-p,"(%d>%d)",pa[idx[t]],pb[idx[t]]);
            for(int a=0;a<m;a++){
                { u64 M=U[a]; if(M && !(M&(M-1))){ long long c0=(long long)Hval; int j=__builtin_ctzll(M); nroots++;
                    if(c0>0&&c0<CMAX&&j<bestj[c0]){ bestj[c0]=j; snprintf(bestdesc[c0],96,"m=%d H=%llu back %s O={%d}",m,Hval,buf,a);} } }
                for(int b=a+1;b<m;b++){
                    u64 w2=W[a][b]+W[b][a]; u64 M=U[a]+U[b];
                    if(!(M&(M-1))){ long long c0=(long long)Hval-(long long)w2; int j=__builtin_ctzll(M); nroots++;
                        if(c0>0&&c0<CMAX&&j<bestj[c0]){ bestj[c0]=j; snprintf(bestdesc[c0],96,"m=%d H=%llu back %s O={%d,%d}",m,Hval,buf,a,b);} }
                    for(int c=b+1;c<m;c++){
                        u64 w3=w2+W[a][c]+W[c][a]+W[b][c]+W[c][b]; u64 M3=M+U[c];
                        if(!(M3&(M3-1))){ long long c0=(long long)Hval-(long long)w3; int j=__builtin_ctzll(M3); nroots++;
                            if(c0>0&&c0<CMAX&&j<bestj[c0]){ bestj[c0]=j; snprintf(bestdesc[c0],96,"m=%d H=%llu back %s O={%d,%d,%d}",m,Hval,buf,a,b,c);} }
                    }
                }
            }
        }
        int t=k-1; while(t>=0 && idx[t]==np-k+t) t--; if(t<0) break; idx[t]++; for(int s=t+1;s<k;s++) idx[s]=idx[s-1]+1;
        if(k==0) break;
    }
    int have=0, full0=0, tot=0; for(int c=1;c<=cmax;c+=2){ tot++; if(bestj[c]<99){ have++; if(bestj[c]==0) full0++; } }
    printf("m=%d back arcs=%d tournaments=%ld cone roots=%ld : odd offsets r <= %d with a cone r + 2^(j+t): %d of %d\n",m,k,nB,nroots,cmax,have,tot);
    
    for(int c=1;c<=cmax;c+=2) if(bestj[c]<99 && (argc>4 || ((c-1)&(c-2))==0)) printf("r=%d  least multiplier 2^%d  %s\n",c,bestj[c],bestdesc[c]);
    return 0;
}
