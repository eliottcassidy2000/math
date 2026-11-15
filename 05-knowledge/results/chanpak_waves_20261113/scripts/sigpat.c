// For each seed with both + and - unit children: multiset of child signatures (P = plus only, M = minus only, B = both, N = none)
#include "armed_lib.h"
static int sig(int n,const int*out){ i64 Pi[MAXN*MAXN],H; build(n,out,Pi,&H); static int p[1<<17],m[1<<17]; int a,b; enum_units(n,Pi,H,p,&a,m,&b,1<<17); return (a>0)+2*(b>0); }
int main(int argc,char**argv){ int n=atoi(argv[1]); CAPN=n+2; size_t NN=(size_t)1<<CAPN;
  he=malloc(sizeof(i64)*NN*CAPN); g=malloc(sizeof(i64)*NN*CAPN); Hc=malloc(sizeof(i64)*NN);
  char line[4096]; const char*S="NPMB";
  static long long pat[4][4][4][4]; // counts of (plus-child has P?,M?,B?) ... simpler: record set masks
  long long table[16][16]={{0}};
  while(fgets(line,sizeof line,stdin)){ if(line[0]=='>') continue;
    int out[MAXN]; parse(line,n,out); i64 Pi[MAXN*MAXN],H; build(n,out,Pi,&H);
    static int p[1<<16],m[1<<16]; int a,b; enum_units(n,Pi,H,p,&a,m,&b,1<<16);
    if(!a||!b||H==1) continue;
    int setp=0,setm=0; int out2[MAXN];
    for(int i=0;i<a;i++){ add_vertex(n,out,p[i],out2); setp|=1<<sig(n+1,out2); }
    for(int i=0;i<b;i++){ add_vertex(n,out,m[i],out2); setm|=1<<sig(n+1,out2); }
    table[setp][setm]++;
  }
  printf("order %d armed seeds (H>1): rows = set of +child signatures, cols = set of -child signatures\n",n);
  for(int i=0;i<16;i++) for(int j=0;j<16;j++) if(table[i][j]){ printf("  +children {"); for(int k=0;k<4;k++) if((i>>k)&1) printf("%c",S[k]); printf("}  -children {"); for(int k=0;k<4;k++) if((j>>k)&1) printf("%c",S[k]); printf("}: %lld\n",table[i][j]); }
  return 0;}
