// betarank.c - the mod-4 law of one-vertex extensions.
// For each tournament R (gentourng format) on m vertices:
//   mu(i,j) = M[i][j] + M[j][i], M[j][i] = #ordered two-path covers (P1 ending at j, P2 starting at i), both non-empty
//   claim: mu(i,j) is even and >= 2;  X -> (H(R+v_X) - H(R))/2 mod 2 is a quadratic form over F_2 with polar form
//   beta(i,j) = mu(i,j)/2 mod 2, so  |#even - #odd| over the 2^m out-sets X is 0 or 2^(m - rank(beta)/2).
// Output: histogram over classes of (rank, exact |bias| numerator), plus consistency failures (must be 0).
#include <stdio.h>
#include <stdlib.h>
typedef unsigned long long u64; typedef long long i64;
#define MAXN 11
static int n; static u64 f[1<<MAXN][MAXN], g[1<<MAXN][MAXN], hS[1<<MAXN]; static int outm[MAXN], inm[MAXN];
int main(int argc,char**argv){ n=atoi(argv[1]); int full=(1<<n)-1; char line[256];
  static u64 M[MAXN][MAXN]; u64 hist[MAXN+1][2]; for(int r=0;r<=MAXN;r++)hist[r][0]=hist[r][1]=0;
  u64 fail_even=0, fail_pred=0, nt=0; u64 minmu=~0ULL; double sumbias=0; i64 signedbias=0;
  while(fgets(line,sizeof line,stdin)){ if(line[0]!='0'&&line[0]!='1')continue; int k=0;
    for(int i=0;i<n;i++){outm[i]=0;inm[i]=0;}
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){ if(line[k++]=='1'){outm[i]|=1<<j;inm[j]|=1<<i;} else {outm[j]|=1<<i;inm[i]|=1<<j;} }
    nt++; hS[0]=1;
    for(int S=1;S<=full;S++){ u64 tot=0; for(int v=0;v<n;v++){ if(!(S>>v&1)){f[S][v]=g[S][v]=0;continue;} int R=S&~(1<<v);
        if(R==0){f[S][v]=g[S][v]=1;tot+=1;continue;} u64 a=0,b=0; int mm=R&inm[v]; while(mm){int w=__builtin_ctz(mm);mm&=mm-1;a+=f[R][w];}
        mm=R&outm[v]; while(mm){int w=__builtin_ctz(mm);mm&=mm-1;b+=g[R][w];} f[S][v]=a;g[S][v]=b;tot+=a; } hS[S]=tot; }
    u64 H=hS[full];
    for(int j=0;j<n;j++)for(int i=0;i<n;i++){ if(i==j){M[j][i]=0;continue;} u64 a=0; int rest=full&~(1<<i)&~(1<<j);
      for(int S0=rest;;S0=(S0-1)&rest){int S=S0|(1<<j); a+=f[S][j]*g[full^S][i]; if(S0==0)break;} M[j][i]=a; }
    int rows[MAXN];
    for(int i=0;i<n;i++){ rows[i]=0; for(int j=0;j<n;j++){ if(i==j)continue; u64 mu=M[i][j]+M[j][i]; if(mu&1)fail_even++; if(mu<minmu)minmu=mu; if((mu>>1)&1) rows[i]|=1<<j; } }
    int rk=0; int rr[MAXN]; for(int i=0;i<n;i++)rr[i]=rows[i];
    for(int c=0;c<n;c++){ int p=-1; for(int i=rk;i<n;i++) if(rr[i]>>c&1){p=i;break;} if(p<0)continue; int t=rr[p];rr[p]=rr[rk];rr[rk]=t;
      for(int i=0;i<n;i++) if(i!=rk && (rr[i]>>c&1)) rr[i]^=rr[rk]; rk++; }
    i64 diff=0;
    for(int X=0;X<=full;X++){ u64 val=0; for(int i=0;i<n;i++){ if(X>>i&1) val+=g[full][i]; else { val+=f[full][i]; int mm=X; while(mm){int t=__builtin_ctz(mm);mm&=mm-1;val+=M[i][t];} } }
      if(((val-H)>>1)&1) diff--; else diff++; }
    i64 ad = diff<0?-diff:diff; i64 pred = 1LL<<(n-rk/2);
    if(!(ad==0 || ad==pred) || (rk&1)) fail_pred++;
    hist[rk][ad!=0]++; sumbias += (double)ad/(double)(1LL<<n); signedbias += diff;
  }
  printf("order %d: classes %llu, parity failures %llu, prediction failures %llu, min mu %llu\n", n, nt, fail_even, fail_pred, minmu);
  for(int r=0;r<=n;r++) if(hist[r][0]||hist[r][1]) printf("  rank %d: %llu classes with zero bias, %llu with bias 2^-%d\n", r, hist[r][0], hist[r][1], r/2);
  printf("  mean |bias| over classes = %.6f ; signed mean (even minus odd)/2^m = %.6f\n", sumbias/nt, (double)signedbias/((double)nt*(double)(1LL<<n)));
  return 0; }
