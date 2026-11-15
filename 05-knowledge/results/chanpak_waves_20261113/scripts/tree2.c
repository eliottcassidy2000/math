// Adaptive full binary tree of unit extensions; a branch whose target is 7 or 21 is not required.
#include "armed_lib.h"
static int armedX(int n,const int*out,int d){
  if(d==0) return 1;
  if(n+1>CAPN) return 0;
  if(node_budget && nodes_visited>=node_budget) return 0;
  nodes_visited++;
  i64 Pi[MAXN*MAXN],H; build(n,out,Pi,&H);
  int *plus=malloc(sizeof(int)*65536),*minus=malloc(sizeof(int)*65536),np,nm;
  enum_units(n,Pi,H,plus,&np,minus,&nm,65536);
  int needp = !(2*H+1==7||2*H+1==21), needm = !(2*H-1==7||2*H-1==21);
  int ok=1,out2[MAXN];
  if(needp){ ok=0; for(int i=0;i<np&&!ok;i++){ add_vertex(n,out,plus[i],out2); if(armedX(n+1,out2,d-1)) ok=1; } }
  if(ok&&needm){ ok=0; for(int i=0;i<nm&&!ok;i++){ add_vertex(n,out,minus[i],out2); if(armedX(n+1,out2,d-1)) ok=1; } }
  free(plus); free(minus); return ok; }
int main(int argc,char**argv){
  int n=atoi(argv[1]),D=atoi(argv[2]); CAPN=n+D+1; if(CAPN>17) CAPN=17;
  if(getenv("BUDGET")) node_budget=atoll(getenv("BUDGET"));
  size_t NN=(size_t)1<<CAPN; he=malloc(sizeof(i64)*NN*CAPN); g=malloc(sizeof(i64)*NN*CAPN); Hc=malloc(sizeof(i64)*NN);
  char line[4096]; int idx=0;
  while(fgets(line,sizeof line,stdin)){
    if(line[0]=='>') continue; idx++;
    int out[MAXN]; parse(line,n,out); i64 H=direct_H(n,out); if(H==1) continue;
    int best=0; long long nv=0;
    for(int d=1;d<=D;d++){ nodes_visited=0; if(armedX(n,out,d)) {best=d; nv=nodes_visited;} else break; }
    line[strcspn(line,"\n")]=0;
    printf("class %d H=%lld depth=%d nodes=%lld code=%s\n",idx,H,best,nv,line); fflush(stdout); }
  return 0;}
