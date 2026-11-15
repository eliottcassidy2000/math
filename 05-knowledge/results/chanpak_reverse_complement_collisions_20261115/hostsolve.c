// Exact host solver for prescribed absolute counts: find a strong marked core C and tapes X, Y >= 1 with
//   F_C(X,Y) = qXY + tX + sY + a(p(X)+p(Y)) + n + f[X odd] + e[Y odd] = N          (THM-5192 (A))
// Each hit is re-derived by the independent cut-transfer recurrence (guard, X-tape, block, Y-tape, guard).
// usage: hostsolve formsfile maxhits_per_target window_cap  < targets
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef __int128 i128; typedef long long i64; typedef unsigned long long u64;
typedef struct { i64 a,q,s,e,t,f,n; int ord; char core[160]; } form;
static i128 p(i128 u){ i128 r=0; while(u>1){ i128 w=u>>1; if(!(u&1)) r+=w+(w&1); u=w; } return r; }
static void pr(i128 x){ char b[64]; int k=0; if(x==0){printf("0");return;} int neg=x<0; if(neg)x=-x; while(x){b[k++]='0'+(int)(x%10);x/=10;} if(neg)putchar('-'); while(k)putchar(b[--k]); }
static i128 rd(const char*s){ i128 x=0; while(*s>='0'&&*s<='9'){x=x*10+(*s-'0');s++;} return x; }
// cut transfer state (h,N,S,d); singleton with mark mu; block C
typedef struct { i128 h,N,S,d; } st;
static st step1(st x,int mu){ st y; y.h=x.h; y.N=x.N+(1-mu)*(x.S+x.d); y.S=2*x.S+mu*x.h; y.d=mu*x.h; return y; }
static st stepb(st x,form*C){ st y; y.h=C->a*x.h; y.N=C->a*x.N+(C->n-C->a)*x.h+C->t*x.S+C->f*x.d; y.S=C->q*x.S+x.h*C->s; y.d=x.h*C->e; return y; }
static i128 transfer(form*C,i128 X,i128 Y){
  st s={1,1,0,0}; s=step1(s,1);                      // marked guard
  int L=0; for(i128 t=X;t>1;t>>=1)L++;
  for(int i=L-1;i>=0;i--) s=step1(s,(int)((X>>i)&1)); // X-tape, MSB first, mark = digit
  s=stepb(s,C);
  int K=0; for(i128 t=Y;t>1;t>>=1)K++;
  for(int j=1;j<=K;j++) s=step1(s,((Y>>(j-1))&1)?0:1); // Y-tape: j-th unmarked iff bit j-1 of Y is 1
  s=step1(s,0);                                       // unmarked guard
  return s.N; }
int main(int argc,char**argv){
  FILE*fp=fopen(argv[1],"r"); int maxhits=atoi(argv[2]); i64 wcap=atoll(argv[3]);
  form F[400]; int nf=0; char line[512];
  while(fgets(line,sizeof line,fp)){ form c; if(sscanf(line,"%lld %lld %lld %lld %lld %lld %lld | %d",&c.a,&c.q,&c.s,&c.e,&c.t,&c.f,&c.n,&c.ord)<8) continue;
    char*bar=strchr(line,'|'); strncpy(c.core,bar+2,159); c.core[strcspn(c.core,"\n")]=0;
    int dup=0; for(int i=0;i<nf;i++) if(F[i].a==c.a&&F[i].q==c.q&&F[i].s==c.s&&F[i].e==c.e&&F[i].t==c.t&&F[i].f==c.f&&F[i].n==c.n) dup=1;
    if(!dup) F[nf++]=c; }
  fprintf(stderr,"%d distinct forms\n",nf);
  char tl[128];
  while(fgets(tl,sizeof tl,stdin)){
    i128 N=rd(tl); if(N<=0) continue;
    int hits=0; long long work=0;
    for(int i=0;i<nf && hits<maxhits;i++){
      form*C=&F[i];
      double Nd=(double)N; double X0=floor(sqrt((double)C->a*Nd/(double)wcap)/(double)C->q); if(X0<1)X0=1;
      double X1=4.0*sqrt(Nd/(double)C->q)+2;
      for(i128 X=(i128)X0; X<=(i128)X1 && hits<maxhits; X++){
        i128 c=C->q*X+C->s; i128 M=N-C->t*X-C->a*p(X)-C->n-C->f*(X&1);
        if(M<c) break;                              // need Y>=1: M >= c*1 + ...
        i128 ylo=(M-C->e)/(c+C->a); if(ylo<1) ylo=1; i128 yhi=M/c;
        for(i128 Y=ylo; Y<=yhi; Y++){ work++;
          if(c*Y+C->a*p(Y)+C->e*(Y&1)==M){
            i128 v=C->q*X*Y+C->t*X+C->s*Y+C->a*(p(X)+p(Y))+C->n+C->f*(X&1)+C->e*(Y&1);
            i128 w=transfer(C,X,Y);
            if(v!=N||w!=N){ printf("VERIFY FAIL\n"); return 1; }
            printf("N="); pr(N); printf(" core(order %d; n,V,out)=[%s] X=",C->ord,C->core); pr(X); printf(" Y="); pr(Y);
            int Lx=0; for(i128 t=X;t>1;t>>=1)Lx++; int Ly=0; for(i128 t=Y;t>1;t>>=1)Ly++;
            printf(" vertices=%d\n",2+Lx+C->ord+Ly+1);
            hits++; if(hits>=maxhits) break; } } } }
    if(!hits){ printf("N="); pr(N); printf(" NO HIT in scanned region (work %lld)\n",work); }
    fflush(stdout);
  }
  return 0; }
