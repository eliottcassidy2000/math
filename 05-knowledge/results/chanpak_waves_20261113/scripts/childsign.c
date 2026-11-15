// For each class (gentourng order n) with both unit extensions, examine every +child and -child:
// does the child itself have a +unit extension, a -unit extension?
#include "armed_lib.h"
int main(int argc,char**argv){
  int n=atoi(argv[1]); CAPN=n+2; size_t NN=(size_t)1<<CAPN;
  he=malloc(sizeof(i64)*NN*CAPN); g=malloc(sizeof(i64)*NN*CAPN); Hc=malloc(sizeof(i64)*NN);
  char line[4096]; long long cnt[2][4]={{0}}; // [sign of move][child type: 0 none,1 plus only,2 minus only,3 both]
  long long seeds=0;
  while(fgets(line,sizeof line,stdin)){
    if(line[0]=='>') continue;
    int out[MAXN]; parse(line,n,out); i64 Pi[MAXN*MAXN],H; build(n,out,Pi,&H);
    if(H==3||H==11) continue;
    static int plus[1<<16],minus[1<<16]; int np,nm; enum_units(n,Pi,H,plus,&np,minus,&nm,1<<16);
    if(!np||!nm) continue; seeds++;
    for(int s=0;s<2;s++){ int *L= s?minus:plus; int m= s?nm:np;
      for(int i=0;i<m;i++){ int out2[MAXN]; add_vertex(n,out,L[i],out2); i64 Pi2[MAXN*MAXN],H2; build(n+1,out2,Pi2,&H2);
        static int p2[1<<16],m2[1<<16]; int a,b; enum_units(n+1,Pi2,H2,p2,&a,m2,&b,1<<16);
        cnt[s][(a>0)+2*(b>0)]++; } }
  }
  printf("order %d armed seeds %lld\n",n,seeds);
  const char*nm[4]={"none","plus-only","minus-only","both"};
  for(int s=0;s<2;s++){ printf(" %s-children:",s?"minus":"plus"); for(int t=0;t<4;t++) printf(" %s=%lld",nm[t],cnt[s][t]); printf("\n"); }
  return 0;}
