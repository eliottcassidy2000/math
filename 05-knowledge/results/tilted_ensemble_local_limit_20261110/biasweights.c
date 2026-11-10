// biasweights.c - exact mod-4 bias of every vertex cube with labelled weights (README section 4.1, deletion-pair form).
// stdin: gentourng output for order m (upper triangle, '1' = i->j for i<j).  Usage: biasweights m
// For each class R on m vertices:
//   f, g, M as in betarank.c;  mu(i,j) = M[i][j] + M[j][i];  beta = mu/2 mod 2;
//   a_i = H(R + v_{i}) - H(R) = g(i) - f(i) + sum_{j != i} M[j][i];  ell_i = a_i/2 mod 2;
//   q(X) = sum_{i in X} ell_i + sum_{{i,j} in X} beta_ij  (Theorem C (2));
//   diff_R = sum_X (-1)^q(X) over all 2^m out-sets, by a Gray code (no use of Theorem C (3));
//   rank(beta) and a basis of its radical by elimination; Theorem C (3) is CHECKED against diff_R;
//   |Aut(R)| by backtracking; labelled weight w_R = m!/|Aut(R)|.
// Output: one line of totals (parsed by merge_bias.py).
#include <stdio.h>
#include <stdlib.h>
typedef unsigned long long u64; typedef long long i64;
#define MAXN 11
static int n; static u64 f[1<<MAXN][MAXN], g[1<<MAXN][MAXN]; static int outm[MAXN], inm[MAXN], sc[MAXN];
static int perm[MAXN], used; static u64 autcount;
static void autrec(int k){
  if(k==n){ autcount++; return; }
  for(int c=0;c<n;c++){ if(used>>c&1) continue; if(sc[c]!=sc[k]) continue; int ok=1;
    for(int j=0;j<k;j++){ int a=(outm[k]>>j)&1, b=(outm[c]>>perm[j])&1; if(a!=b){ok=0;break;} }
    if(!ok) continue; perm[k]=c; used|=1<<c; autrec(k+1); used&=~(1<<c); }
}
int main(int argc,char**argv){ n=atoi(argv[1]); int full=(1<<n)-1; char line[256]; static u64 M[MAXN][MAXN];
  u64 fact=1; for(int i=2;i<=n;i++) fact*=i;
  u64 cls[MAXN+1][2], lab[MAXN+1][2]; for(int r=0;r<=MAXN;r++){cls[r][0]=cls[r][1]=lab[r][0]=lab[r][1]=0;}
  u64 nt=0, wtot=0, fail_even=0, fail_a=0, fail_c3=0, fail_aut=0, absum=0; i64 signedS=0, plainU=0;
  while(fgets(line,sizeof line,stdin)){ if(line[0]!='0'&&line[0]!='1')continue; int k=0;
    for(int i=0;i<n;i++){outm[i]=0;inm[i]=0;}
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){ if(line[k++]=='1'){outm[i]|=1<<j;inm[j]|=1<<i;} else {outm[j]|=1<<i;inm[i]|=1<<j;} }
    nt++;
    for(int S=1;S<=full;S++){ for(int v=0;v<n;v++){ if(!(S>>v&1)){f[S][v]=g[S][v]=0;continue;} int R=S&~(1<<v);
        if(R==0){f[S][v]=g[S][v]=1;continue;} u64 a=0,b=0; int mm=R&inm[v]; while(mm){int w=__builtin_ctz(mm);mm&=mm-1;a+=f[R][w];}
        mm=R&outm[v]; while(mm){int w=__builtin_ctz(mm);mm&=mm-1;b+=g[R][w];} f[S][v]=a;g[S][v]=b; } }
    u64 H=0; for(int v=0;v<n;v++) H+=f[full][v];
    for(int j=0;j<n;j++)for(int i=0;i<n;i++){ if(i==j){M[j][i]=0;continue;} u64 a=0; int rest=full&~(1<<i)&~(1<<j);
      for(int S0=rest;;S0=(S0-1)&rest){int S=S0|(1<<j); a+=f[S][j]*g[full^S][i]; if(S0==0)break;} M[j][i]=a; }
    int rows[MAXN], ell[MAXN];
    for(int i=0;i<n;i++){ rows[i]=0; for(int j=0;j<n;j++){ if(i==j)continue; u64 mu=M[i][j]+M[j][i]; if(mu&1)fail_even++; if((mu>>1)&1) rows[i]|=1<<j; }
      i64 a=(i64)g[full][i]-(i64)f[full][i]; for(int j=0;j<n;j++) if(j!=i) a+=(i64)M[j][i]; if(a&1) fail_a++; ell[i]=(int)((a>>1)&1); }
    // rank and radical: eliminate on [beta | I]
    int rr[MAXN], id[MAXN]; for(int i=0;i<n;i++){rr[i]=rows[i]; id[i]=1<<i;} int rk=0;
    for(int c=0;c<n;c++){ int p=-1; for(int i=rk;i<n;i++) if(rr[i]>>c&1){p=i;break;} if(p<0)continue;
      int t=rr[p];rr[p]=rr[rk];rr[rk]=t; t=id[p];id[p]=id[rk];id[rk]=t;
      for(int i=0;i<n;i++) if(i!=rk && (rr[i]>>c&1)){ rr[i]^=rr[rk]; id[i]^=id[rk]; } rk++; }
    int radzero=1;   // does q vanish on the radical?  (q is additive there: test a basis)
    for(int i=rk;i<n;i++){ int X=id[i]; int q=0; for(int a=0;a<n;a++) if(X>>a&1){ q^=ell[a]; q^=__builtin_popcount(rows[a]&X&((1<<a)-1))&1; } if(q){radzero=0;} }
    // exact Gauss sum by Gray code
    i64 diff=1; { int X=0, q=0; for(int t=1;t<=full;t++){ int i=__builtin_ctz(t); q^=ell[i]^(__builtin_popcount(rows[i]&X&~(1<<i))&1); X^=1<<i; diff += q?-1:1; } }
    i64 ad=diff<0?-diff:diff; i64 pred=1LL<<(n-rk/2);
    if((rk&1) || !((ad==0&&!radzero)||(ad==pred&&radzero))) fail_c3++;
    for(int i=0;i<n;i++) sc[i]=__builtin_popcount(outm[i]); autcount=0; used=0; autrec(0);
    if(autcount==0 || fact%autcount) { fail_aut++; continue; }
    u64 w=fact/autcount; wtot+=w; cls[rk][ad!=0]++; lab[rk][ad!=0]+=w;
    i64 eps=(((H-1)>>1)&1)?-1:1; signedS+=eps*(i64)w*diff; plainU+=(i64)w*diff; absum+=w*(u64)ad;
  }
  printf("order %d classes %llu weight %llu failures(mu odd, a odd, C3, aut) %llu %llu %llu %llu S %lld U %lld A %llu", n, nt, wtot, fail_even, fail_a, fail_c3, fail_aut, signedS, plainU, absum);
  for(int r=0;r<=n;r++) printf(" r%d %llu %llu %llu %llu", r, cls[r][0], cls[r][1], lab[r][0], lab[r][1]);
  printf("\n"); return 0; }
