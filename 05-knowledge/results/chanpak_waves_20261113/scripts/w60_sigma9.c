/* Smoothness sigma (largest prime-quotient order) of every order-9 tournament (one-vertex extensions of the 6,880
   order-8 classes), together with orders <= 8 from the class list: minimal sigma per count h <= LIM. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static u64 Hcount(int m, const int*o){ static u64 dp[512][9]; for(int S=0;S<(1<<m);S++) for(int v=0;v<m;v++) dp[S][v]=0; for(int v=0;v<m;v++) dp[1<<v][v]=1;
  for(int S=1;S<(1<<m);S++) for(int v=0;v<m;v++) if(dp[S][v]){ int t=o[v]&~S; while(t){ int w=__builtin_ctz(t); t&=t-1; dp[S|1<<w][w]+=dp[S][v]; } }
  u64 h=0; for(int v=0;v<m;v++) h+=dp[(1<<m)-1][v]; return h; }
static int ismodule(int M,int V,const int*o){ for(int x=0;x<9;x++) if((V>>x&1)&&!(M>>x&1)){ int a=-1; for(int y=0;y<9;y++) if(M>>y&1){ int b=o[x]>>y&1; if(a<0) a=b; else if(a!=b) return 0; } } return 1; }
static int reachset(int s,int V,const int*o,int fwd){ int r=1<<s,fr=r; while(fr){ int nf=0; for(int v=0;v<9;v++) if(fr>>v&1){ if(fwd) nf|=o[v]&V; else for(int u=0;u<9;u++) if((V>>u&1)&&(o[u]>>v&1)) nf|=1<<u; } nf&=~r; r|=nf; fr=nf; } return r; }
static int sigma(int V,const int*o){ int c=__builtin_popcount(V); if(c<=1) return 1;
  int s=__builtin_ctz(V); if((reachset(s,V,o,1)&reachset(s,V,o,0))!=V || reachset(s,V,o,1)!=V){ /* not strong: recurse on SCCs */
    int rem=V,best=1; while(rem){ int t=__builtin_ctz(rem); int scc=reachset(t,V,o,1)&reachset(t,V,o,0); int g=sigma(scc,o); if(g>best) best=g; rem&=~scc; } return best; }
  static int mods[9][512]; int depth=0; (void)depth; int ml[512], nm=0;
  for(int M=(V-1)&V; M; M=(M-1)&V){ if(ismodule(M,V,o)) ml[nm++]=M; }
  int parts[9],np=0,cover=0;
  for(int i=0;i<nm;i++){ int M=ml[i]; int maximal=1; for(int j=0;j<nm;j++) if(ml[j]!=M&&(ml[j]&M)==M){ maximal=0; break; } if(maximal&&!(cover&M)){ parts[np++]=M; cover|=M; } }
  for(int v=0;v<9;v++) if((V>>v&1)&&!(cover>>v&1)){ parts[np++]=1<<v; cover|=1<<v; }
  int best=np; for(int i=0;i<np;i++){ int g=sigma(parts[i],o); if(g>best) best=g; } return best; }
int main(int argc,char**argv){ int LIM=atoi(argv[1]); int *minsig=malloc(sizeof(int)*(LIM+1)); for(int h=0;h<=LIM;h++) minsig[h]=99;
  FILE*fi=fopen("w59_reps8.txt","r"); char line[256]; static int reps8[7000][8]; int k8=0;
  while(fgets(line,sizeof line,fi)){ char*p=line; int n=strtol(p,&p,10); int o[9]; for(int i=0;i<n;i++) o[i]=strtol(p,&p,10);
    u64 h=Hcount(n,o); int g=sigma((1<<n)-1,o); if(h<=(u64)LIM && g<minsig[h]) minsig[h]=g; if(n==8) memcpy(reps8[k8++],o,sizeof(int)*8); }
  fclose(fi);
  for(int r=0;r<k8;r++) for(int S=0;S<256;S++){ int o[9]; for(int i=0;i<8;i++){ o[i]=reps8[r][i]; if(!(S>>i&1)) o[i]|=1<<8; } o[8]=S;
      u64 h=Hcount(9,o); if(h>(u64)LIM) continue; int g=sigma(511,o); if(g<minsig[h]) minsig[h]=g; }
  for(int h=1;h<=LIM;h+=2) if(minsig[h]<99) printf("%d %d\n",h,minsig[h]);
  return 0; }
