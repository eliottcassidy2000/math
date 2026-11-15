// For armed seeds: number of distinct child types per sign, using an invariant (H, sorted scores, sorted per-vertex
// (Pi row sum, Pi col sum)) that is converse-symmetrized.
#include "armed_lib.h"
typedef unsigned long long u64;
static u64 inv(int n,const int*out){ i64 Pi[MAXN*MAXN],H; build(n,out,Pi,&H);
  i64 key[MAXN]; for(int v=0;v<n;v++){ i64 r=0,c=0; for(int u=0;u<n;u++){ r+=Pi[v*n+u]; c+=Pi[u*n+v]; } int sc=__builtin_popcount(out[v]);
    i64 a=r*1000003+c, b=c*1000003+r; key[v]=(a<b?a:b)*64+ (sc<n-1-sc?sc:n-1-sc); }
  for(int i=0;i<n;i++) for(int j=i+1;j<n;j++) if(key[j]<key[i]){ i64 t=key[i]; key[i]=key[j]; key[j]=t; }
  u64 h=1469598103934665603ULL^(u64)H; for(int v=0;v<n;v++){ h^=(u64)key[v]; h*=1099511628211ULL; } return h; }
int main(int argc,char**argv){ int n=atoi(argv[1]); CAPN=n+2; size_t NN=(size_t)1<<CAPN;
  he=malloc(sizeof(i64)*NN*CAPN); g=malloc(sizeof(i64)*NN*CAPN); Hc=malloc(sizeof(i64)*NN);
  char line[4096]; long long hist[2][8]={{0}}; long long seeds=0;
  while(fgets(line,sizeof line,stdin)){ if(line[0]=='>') continue;
    int out[MAXN]; parse(line,n,out); i64 Pi[MAXN*MAXN],H; build(n,out,Pi,&H);
    static int p[1<<16],m[1<<16]; int a,b; enum_units(n,Pi,H,p,&a,m,&b,1<<16);
    if(!a||!b||H==1) continue; seeds++;
    for(int s=0;s<2;s++){ int*L=s?m:p; int c=s?b:a; u64 seen[4096]; int ns=0; int out2[MAXN];
      for(int i=0;i<c;i++){ add_vertex(n,out,L[i],out2); u64 k=inv(n+1,out2); int f=0; for(int j=0;j<ns;j++) if(seen[j]==k) f=1; if(!f) seen[ns++]=k; }
      hist[s][ns<7?ns:7]++; }
  }
  printf("order %d armed seeds %lld; distinct +child types histogram:",n,seeds); for(int k=1;k<8;k++) printf(" %d:%lld",k,hist[0][k]);
  printf("; -child types:"); for(int k=1;k<8;k++) printf(" %d:%lld",k,hist[1][k]); printf("\n"); return 0;}
