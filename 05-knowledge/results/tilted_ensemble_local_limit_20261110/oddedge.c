// oddedge.c - Lemma 4.5 of README.md: for every tournament R on m >= 2 vertices
//     sum_{i<j} mu_R(i,j) = sum_i a_i = sum_{S proper, non-empty} H(R[S]) H(R[V \ S]) = 2 (mod 4),
// so beta_R = mu_R/2 mod 2 has an odd number of entries above the diagonal and the linear part ell = a/2 mod 2 has
// odd weight.  stdin: gentourng output for order m.  Usage: oddedge m
#include <stdio.h>
#include <stdlib.h>
typedef unsigned long long u64; typedef long long i64;
#define MAXN 11
static int n; static u64 f[1<<MAXN][MAXN], g[1<<MAXN][MAXN], hS[1<<MAXN]; static int outm[MAXN], inm[MAXN];
int main(int argc,char**argv){ n=atoi(argv[1]); int full=(1<<n)-1; char line[256]; static u64 M[MAXN][MAXN];
  u64 nt=0, f_id1=0, f_id2=0, f_mod4=0, f_edges=0, f_ell=0; u64 edgehist[64]; for(int i=0;i<64;i++)edgehist[i]=0;
  while(fgets(line,sizeof line,stdin)){ if(line[0]!='0'&&line[0]!='1')continue; int k=0;
    for(int i=0;i<n;i++){outm[i]=0;inm[i]=0;}
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){ if(line[k++]=='1'){outm[i]|=1<<j;inm[j]|=1<<i;} else {outm[j]|=1<<i;inm[i]|=1<<j;} }
    nt++; hS[0]=1;
    for(int S=1;S<=full;S++){ u64 tot=0; for(int v=0;v<n;v++){ if(!(S>>v&1)){f[S][v]=g[S][v]=0;continue;} int R=S&~(1<<v);
        if(R==0){f[S][v]=g[S][v]=1;tot+=1;continue;} u64 a=0,b=0; int mm=R&inm[v]; while(mm){int w=__builtin_ctz(mm);mm&=mm-1;a+=f[R][w];}
        mm=R&outm[v]; while(mm){int w=__builtin_ctz(mm);mm&=mm-1;b+=g[R][w];} f[S][v]=a;g[S][v]=b;tot+=a; } hS[S]=tot; }
    for(int j=0;j<n;j++)for(int i=0;i<n;i++){ if(i==j){M[j][i]=0;continue;} u64 a=0; int rest=full&~(1<<i)&~(1<<j);
      for(int S0=rest;;S0=(S0-1)&rest){int S=S0|(1<<j); a+=f[S][j]*g[full^S][i]; if(S0==0)break;} M[j][i]=a; }
    u64 summu=0, suma=0, sumsplit=0; int edges=0, ellw=0;
    for(int i=0;i<n;i++){ i64 a=(i64)g[full][i]-(i64)f[full][i]; for(int j=0;j<n;j++){ if(j==i)continue; a+=(i64)M[j][i];
        if(j>i){ u64 mu=M[i][j]+M[j][i]; summu+=mu; edges+=(int)((mu>>1)&1); } }
      suma+=(u64)a; ellw+=(int)((a>>1)&1); }
    for(int S=1;S<full;S++) sumsplit+=hS[S]*hS[full^S];
    if(summu!=suma) f_id1++; if(suma!=sumsplit) f_id2++; if((summu&3)!=2) f_mod4++; if(!(edges&1)) f_edges++; if(!(ellw&1)) f_ell++;
    edgehist[edges]++;
  }
  printf("order %d classes %llu failures: sum mu != sum a %llu ; sum a != sum_S H(S)H(S^c) %llu ; sum mu != 2 mod 4 %llu ; even number of beta-edges %llu ; even weight of ell %llu ; edges:", n, nt, f_id1, f_id2, f_mod4, f_edges, f_ell);
  for(int e=0;e<64;e++) if(edgehist[e]) printf(" %d:%llu", e, edgehist[e]);
  printf("\n"); return 0; }
