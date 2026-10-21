/* tour_census.c -- per-class invariants of all tournaments of order n (from gentourng ascii).
 * mac-mini-2026-10-21-fanopsp.
 * Input: gentourng default output, one line per class: upper triangle row by row,
 *   character for pair (i,j), i<j, is '1' iff i -> j.
 * Output per class: H (Hamiltonian paths), C (directed Hamiltonian cycles), strong (0/1),
 *   aut (|Aut|), c3 (cyclic triangles), selfconv (0/1).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int n; static uint32_t out[32], in_[32];
static int adj(int i,int j){ return (out[i]>>j)&1; }
static uint64_t f[1<<16][16];
static uint64_t hampaths(void){
  int full=(1<<n)-1; for(int m=0;m<=full;m++) for(int v=0;v<n;v++) f[m][v]=0;
  for(int v=0;v<n;v++) f[1<<v][v]=1;
  for(int m=1;m<=full;m++) for(int v=0;v<n;v++){ uint64_t c=f[m][v]; if(!c) continue;
    uint32_t cand=out[v]&~m; while(cand){int w=__builtin_ctz(cand); cand&=cand-1; f[m|(1<<w)][w]+=c;} }
  uint64_t H=0; for(int v=0;v<n;v++) H+=f[full][v]; return H; }
static uint64_t hamcycles(void){ /* directed Ham cycles through vertex 0 */
  if(n<3) return 0; int full=(1<<n)-1; for(int m=0;m<=full;m++) for(int v=0;v<n;v++) f[m][v]=0;
  f[1][0]=1;
  for(int m=1;m<=full;m+=2) for(int v=0;v<n;v++){ uint64_t c=f[m][v]; if(!c) continue;
    uint32_t cand=out[v]&~m; while(cand){int w=__builtin_ctz(cand); cand&=cand-1; f[m|(1<<w)][w]+=c;} }
  uint64_t C=0; for(int v=1;v<n;v++) if(adj(v,0)) C+=f[full][v]; return C; }
static int strong(void){ uint32_t r=1, prev=0; while(r!=prev){prev=r; for(int v=0;v<n;v++) if(r>>v&1) r|=out[v];}
  if(r!=(1u<<n)-1) return 0; r=1;prev=0; while(r!=prev){prev=r; for(int v=0;v<n;v++) if(r>>v&1) r|=in_[v];}
  return r==(1u<<n)-1; }
/* automorphisms (anti=0) / anti-automorphisms (anti=1) by backtracking with invariant colours */
static int col[32], tgt[32], img[32], used[32]; static long cnt; static int anti;
static int ok(int v,int w){ for(int u=0;u<v;u++){ int a=adj(u,v), b=adj(img[u],w); if(anti) b=!b; if(a!=b) return 0; } return 1; }
static void bt(int v){ if(v==n){cnt++;return;}
  for(int w=0;w<n;w++) if(!used[w] && col[w]==tgt[v] && ok(v,w)){ used[w]=1; img[v]=w; bt(v+1); used[w]=0; } }
static long hcol(int s,const int*ho,const int*hi){ long h=s; for(int k=0;k<n;k++) h=(h*131+ho[k])%1000000007L; for(int k=0;k<n;k++) h=(h*137+hi[k])%1000000007L; return h; }
static long countmaps(int a){ anti=a; int sc[32]; for(int v=0;v<n;v++) sc[v]=__builtin_popcount(out[v]);
  for(int v=0;v<n;v++){ int ho[32]={0},hi[32]={0}; for(int w=0;w<n;w++){ if(w==v) continue; if(adj(v,w)) ho[sc[w]]++; else hi[sc[w]]++; }
    col[v]=(int)hcol(sc[v],ho,hi);
    if(!a) tgt[v]=col[v];
    else { int ho2[32]={0},hi2[32]={0}; for(int w=0;w<n;w++){ if(w==v) continue; if(adj(w,v)) ho2[n-1-sc[w]]++; else hi2[n-1-sc[w]]++; }
      tgt[v]=(int)hcol(n-1-sc[v],ho2,hi2); } }
  cnt=0; memset(used,0,sizeof used); bt(0); return cnt; }
int main(int argc,char**argv){ n=atoi(argv[1]); char line[4096];
  while(fgets(line,sizeof line,stdin)){ if(line[0]!='0'&&line[0]!='1') continue;
    memset(out,0,sizeof out); memset(in_,0,sizeof in_); int p=0;
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if(line[p]=='1'){out[i]|=1u<<j; in_[j]|=1u<<i;} else {out[j]|=1u<<i; in_[i]|=1u<<j;} p++; }
    uint64_t H=hampaths(), C=hamcycles(); int s=strong(); long A=countmaps(0); long AA=countmaps(1);
    int c3=0; for(int a=0;a<n;a++) for(int b=a+1;b<n;b++) for(int c=b+1;c<n;c++){ int x=adj(a,b)+adj(b,c)+adj(c,a); if(x==0||x==3) c3++; }
    printf("%llu %llu %d %ld %ld %d\n",(unsigned long long)H,(unsigned long long)C,s,A,AA,c3); }
  return 0; }
