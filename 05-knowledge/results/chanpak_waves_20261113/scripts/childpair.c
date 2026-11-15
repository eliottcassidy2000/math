// For each armed seed: H, number of + children, armed + children, number of - children, armed - children
#include "armed_lib.h"
static int isarmed(int n,const int*out){ i64 Pi[MAXN*MAXN],H; build(n,out,Pi,&H);
  static int p[1<<17],m[1<<17]; int a,b; enum_units(n,Pi,H,p,&a,m,&b,1<<17); return (a>0)+2*(b>0); }
int main(int argc,char**argv){
  int n=atoi(argv[1]); CAPN=n+2; size_t NN=(size_t)1<<CAPN;
  he=malloc(sizeof(i64)*NN*CAPN); g=malloc(sizeof(i64)*NN*CAPN); Hc=malloc(sizeof(i64)*NN);
  char line[4096]; int idx=0;
  while(fgets(line,sizeof line,stdin)){
    if(line[0]=='>') continue; idx++;
    int out[MAXN]; parse(line,n,out); i64 Pi[MAXN*MAXN],H; build(n,out,Pi,&H);
    static int plus[1<<16],minus[1<<16]; int np,nm; enum_units(n,Pi,H,plus,&np,minus,&nm,1<<16);
    if(!np||!nm) continue;
    int ap=0,am=0,pp=0,pm=0,mp=0,mm=0; int out2[MAXN];
    for(int i=0;i<np;i++){ add_vertex(n,out,plus[i],out2); int t=isarmed(n+1,out2); if(t==3) ap++; if(t&1) pp++; if(t&2) pm++; }
    for(int i=0;i<nm;i++){ add_vertex(n,out,minus[i],out2); int t=isarmed(n+1,out2); if(t==3) am++; if(t&1) mp++; if(t&2) mm++; }
    printf("class %d H=%lld  +children %d (armed %d; with+ %d, with- %d)   -children %d (armed %d; with+ %d, with- %d)\n",idx,H,np,ap,pp,pm,nm,am,mp,mm);
  }
  return 0;}
