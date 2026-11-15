// For each class: H, strong?, #plus units, #minus units; for each unit child: strong?
#include "armed_lib.h"
static int strong(int n,const int*out){ int reach=1,fr=1; while(fr){int nf=0; for(int v=0;v<n;v++) if((fr>>v)&1) nf|=out[v]; nf&=~reach; reach|=nf; fr=nf;}
  if(reach!=(1<<n)-1) return 0; int in[MAXN]; for(int v=0;v<n;v++) in[v]=0; for(int v=0;v<n;v++) for(int u=0;u<n;u++) if((out[v]>>u)&1) in[u]|=1<<v;
  reach=1; fr=1; while(fr){int nf=0; for(int v=0;v<n;v++) if((fr>>v)&1) nf|=in[v]; nf&=~reach; reach|=nf; fr=nf;} return reach==(1<<n)-1; }
int main(int argc,char**argv){ int n=atoi(argv[1]); CAPN=n+2; size_t NN=(size_t)1<<CAPN;
  he=malloc(sizeof(i64)*NN*CAPN); g=malloc(sizeof(i64)*NN*CAPN); Hc=malloc(sizeof(i64)*NN);
  char line[4096]; long long cls=0, armed=0, armed_nonstrong=0, units=0, unit_child_strong=0, plus_only_nonstrong=0, minus_only_nonstrong=0, one_sided=0;
  long long redwithunits=0;
  while(fgets(line,sizeof line,stdin)){ if(line[0]=='>') continue; cls++;
    int out[MAXN]; parse(line,n,out); i64 Pi[MAXN*MAXN],H; build(n,out,Pi,&H);
    static int p[1<<16],m[1<<16]; int a,b; enum_units(n,Pi,H,p,&a,m,&b,1<<16);
    if(H==1) continue;
    int st=strong(n,out);
    if(a&&b){ armed++; if(!st) armed_nonstrong++; }
    else if(a||b){ one_sided++; if(!st){ if(a) plus_only_nonstrong++; else minus_only_nonstrong++; } }
    if(!st && (a||b)) redwithunits++;
    int out2[MAXN];
    for(int s=0;s<2;s++){ int*L=s?m:p; int c=s?b:a; for(int i=0;i<c;i++){ add_vertex(n,out,L[i],out2); units++; if(strong(n+1,out2)) unit_child_strong++; } }
  }
  printf("order %d: classes %lld (H>1 considered); armed %lld, armed but not strong %lld; one-sided %lld (non-strong: plus-only %lld, minus-only %lld); unit children %lld, strong %lld\n",
    n,cls,armed,armed_nonstrong,one_sided,plus_only_nonstrong,minus_only_nonstrong,units,unit_child_strong);
  return 0;}
