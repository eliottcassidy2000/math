// Read gentourng lines (upper-triangle bits, 1 means i->j for i<j).
// For each tournament A and root v output: H(A), H(A-v), L(A,v), E=end_v(A), S=start_v(A)
// L(A,v) = #ordered pairs (P,Q) of disjoint paths covering A-v (empty allowed), Q empty or v->first(Q).
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define MAXN 12
typedef unsigned long long u64;
int n; int adj[MAXN][MAXN];
static u64 endc[1<<MAXN][MAXN], startc[1<<MAXN][MAXN];
// paths in subgraph on vertex set given by list
void paths_dp(int m, int *vs, int (*a)[MAXN]){
  // endc[S][x]: # paths covering S (subset of 0..m-1 indices) ending at x ; startc similarly
  int full=1<<m;
  for(int S=0;S<full;S++) for(int x=0;x<m;x++){endc[S][x]=0;startc[S][x]=0;}
  for(int x=0;x<m;x++){endc[1<<x][x]=1;startc[1<<x][x]=1;}
  for(int S=1;S<full;S++) for(int x=0;x<m;x++){
    if(!(S>>x&1)) continue;
    u64 e=endc[S][x], s=startc[S][x];
    if(!e && !s) continue;
    for(int y=0;y<m;y++){ if(S>>y&1) continue;
      if(e && a[vs[x]][vs[y]]) endc[S|1<<y][y]+=e;   // extend at end: x->y
      if(s && a[vs[y]][vs[x]]) startc[S|1<<y][y]+=s; // extend at start: y->x
    }
  }
}
int main(int argc,char**argv){
  char line[256];
  while(fgets(line,sizeof line,stdin)){
    if(line[0]=='>') continue;
    int len=strlen(line); while(len&&(line[len-1]=='\n'||line[len-1]=='\r')) line[--len]=0;
    // determine n from len = n(n-1)/2
    n=0; while(n*(n-1)/2<len) n++;
    int k=0; for(int i=0;i<n;i++) for(int j=0;j<n;j++) adj[i][j]=0;
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if(line[k]=='1') adj[i][j]=1; else adj[j][i]=1; k++; }
    int vs[MAXN]; for(int i=0;i<n;i++) vs[i]=i;
    paths_dp(n,vs,adj);
    u64 H=0; int full=(1<<n)-1; for(int x=0;x<n;x++) H+=endc[full][x];
    u64 Ev[MAXN],Sv[MAXN]; for(int x=0;x<n;x++){Ev[x]=endc[full][x];Sv[x]=startc[full][x];}
    printf("%s %llu",line,H);
    for(int v=0;v<n;v++){
      int m=0; int w[MAXN]; for(int i=0;i<n;i++) if(i!=v) w[m++]=i;
      paths_dp(m,w,adj);
      int fm=(1<<m)-1; u64 h=0; for(int x=0;x<m;x++) h+=endc[fm][x];
      u64 L=0;
      for(int S=0;S<=fm;S++){
        u64 pcount = (S==0)?1:0; if(S){ for(int x=0;x<m;x++) pcount+=endc[S][x]; }
        if(!pcount) continue;
        int C=fm^S; u64 q;
        if(C==0) q=1; else { q=0; for(int y=0;y<m;y++) if((C>>y&1) && adj[v][w[y]]) q+=startc[C][y]; }
        L+=pcount*q;
      }
      printf(" %d:%llu,%llu,%llu,%llu",v,h,L,Ev[v],Sv[v]);
    }
    printf("\n");
  }
  return 0;
}
