/* anatomy.c -- for strong tournaments of order n (gentourng -c) with H in a target set, print
 * H, alpha_1 (directed odd cycles), alpha_2 (disjoint pairs), alpha_3, c3,c5,c7,c9, score multiset.
 * Checks H = 1 + 2a1 + 4a2 + 8a3.  mac-mini-2026-10-21-fanopsp.  Usage: gentourng -c -q n | ./anatomy n h1 h2 ... */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int n; static uint32_t out[16];
static uint64_t f[1<<12][12];
static uint64_t hampaths(void){ int full=(1<<n)-1; memset(f,0,sizeof(uint64_t)*(full+1)*12);
  for(int v=0;v<n;v++) f[1<<v][v]=1;
  for(int m=1;m<=full;m++) for(int v=0;v<n;v++){ uint64_t c=f[m][v]; if(!c) continue; uint32_t cand=out[v]&~m;
      while(cand){int w=__builtin_ctz(cand); cand&=cand-1; f[m|(1<<w)][w]+=c;} }
  uint64_t H=0; for(int v=0;v<n;v++) H+=f[full][v]; return H; }
static uint64_t cyc[1<<12]; /* directed Ham cycles of T[mask] */
static void allcycles(void){ int full=(1<<n)-1;
  for(int mask=1;mask<=full;mask++){ cyc[mask]=0; int k=__builtin_popcount(mask); if(k<3) continue;
    int s=__builtin_ctz(mask); /* cycles through lowest vertex s, counted once */
    static uint64_t g[1<<12][12]; /* DP restricted to subsets of mask containing s */
    /* simple DP over submasks: iterate over submasks of mask containing s */
    int rest=mask&~(1<<s); /* enumerate subsets of rest */
    /* use a local map from submask index to dp: compress bits */
    int idx[12],ni=0; for(int v=0;v<n;v++) if(rest>>v&1) idx[ni++]=v;
    int lim=1<<ni; static uint64_t h[1<<11][12];
    for(int a=0;a<lim;a++) for(int t=0;t<ni;t++) h[a][t]=0;
    for(int t=0;t<ni;t++) if(out[s]>>idx[t]&1) h[1<<t][t]=1;
    for(int a=1;a<lim;a++) for(int t=0;t<ni;t++){ uint64_t c=h[a][t]; if(!c) continue;
      for(int u=0;u<ni;u++){ if(a>>u&1) continue; if(out[idx[t]]>>idx[u]&1) h[a|(1<<u)][u]+=c; } }
    uint64_t C=0; for(int t=0;t<ni;t++) if(out[idx[t]]>>s&1) C+=h[lim-1][t]; cyc[mask]=C; (void)g; } }
int main(int argc,char**argv){ n=atoi(argv[1]); int nt=argc-2; long tg[64]; for(int i=0;i<nt;i++) tg[i]=atol(argv[i+2]);
  char line[512];
  while(fgets(line,sizeof line,stdin)){ if(line[0]!='0'&&line[0]!='1') continue; memset(out,0,sizeof out); int p=0;
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if(line[p]=='1') out[i]|=1u<<j; else out[j]|=1u<<i; p++; }
    uint64_t H=hampaths(); int hit=0; for(int i=0;i<nt;i++) if((long)H==tg[i]) hit=1; if(!hit) continue;
    allcycles(); int full=(1<<n)-1; uint64_t ck[16]={0}, a1=0,a2=0,a3=0;
    for(int m=1;m<=full;m++){ int k=__builtin_popcount(m); if(k%2==1 && k>=3){ ck[k]+=cyc[m]; a1+=cyc[m]; } }
    for(int m1=1;m1<=full;m1++){ int k1=__builtin_popcount(m1); if(k1%2==0||k1<3||!cyc[m1]) continue;
      int rem=full&~m1; for(int m2=rem;m2;m2=(m2-1)&rem){ if(m2<m1) continue; int k2=__builtin_popcount(m2); if(k2%2==0||k2<3||!cyc[m2]) continue;
        a2+=cyc[m1]*cyc[m2];
        int rem2=rem&~m2; for(int m3=rem2;m3;m3=(m3-1)&rem2){ if(m3<m2) continue; int k3=__builtin_popcount(m3); if(k3%2==0||k3<3||!cyc[m3]) continue; a3+=cyc[m1]*cyc[m2]*cyc[m3]; } } }
    int sc[16]; for(int v=0;v<n;v++) sc[v]=__builtin_popcount(out[v]);
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++) if(sc[j]<sc[i]){int t=sc[i];sc[i]=sc[j];sc[j]=t;}
    printf("H=%llu a1=%llu a2=%llu a3=%llu check=%s c3=%llu c5=%llu c7=%llu c9=%llu scores=",(unsigned long long)H,(unsigned long long)a1,(unsigned long long)a2,(unsigned long long)a3,(1+2*a1+4*a2+8*a3==H)?"ok":"BAD",(unsigned long long)ck[3],(unsigned long long)ck[5],(unsigned long long)ck[7],(unsigned long long)ck[9]);
    for(int v=0;v<n;v++) printf("%d",sc[v]); printf("\n"); }
  return 0; }
