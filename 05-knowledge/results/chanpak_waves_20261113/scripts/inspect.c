// For armed seeds of order n with H>=9: print score sequence-like info and unit out-sets; for each child, its unit out-sets
#include "armed_lib.h"
static void show(int n,const int*out,const char*tag){ i64 Pi[MAXN*MAXN],H; build(n,out,Pi,&H);
  static int p[1<<16],m[1<<16]; int a,b; enum_units(n,Pi,H,p,&a,m,&b,1<<16);
  printf("%s order %d H=%lld outdeg:",tag,n,H); for(int v=0;v<n;v++) printf(" %d",__builtin_popcount(out[v])); printf(" | +X:");
  for(int i=0;i<a;i++){ printf(" {"); for(int v=0;v<n;v++) if((p[i]>>v)&1) printf("%d",v); printf("}"); }
  printf(" -X:"); for(int i=0;i<b;i++){ printf(" {"); for(int v=0;v<n;v++) if((m[i]>>v)&1) printf("%d",v); printf("}"); } printf("\n"); }
int main(int argc,char**argv){ int n=atoi(argv[1]); int want=atoi(argv[2]); CAPN=n+2; size_t NN=(size_t)1<<CAPN;
  he=malloc(sizeof(i64)*NN*CAPN); g=malloc(sizeof(i64)*NN*CAPN); Hc=malloc(sizeof(i64)*NN);
  char line[4096]; int idx=0;
  while(fgets(line,sizeof line,stdin)){ if(line[0]=='>') continue; idx++; if(idx!=want) continue;
    int out[MAXN]; parse(line,n,out); show(n,out,"SEED");
    i64 Pi[MAXN*MAXN],H; build(n,out,Pi,&H); static int p[1<<16],m[1<<16]; int a,b; enum_units(n,Pi,H,p,&a,m,&b,1<<16);
    int out2[MAXN];
    for(int i=0;i<a;i++){ add_vertex(n,out,p[i],out2); char t[64]; sprintf(t,"  +child X=%d",p[i]); show(n+1,out2,t); }
    for(int i=0;i<b;i++){ add_vertex(n,out,m[i],out2); char t[64]; sprintf(t,"  -child X=%d",m[i]); show(n+1,out2,t); }
  }
  return 0;}
