// k35cube.c - the triangle and pentagon parity races on a vertex cube, exactly, by linear algebra over F_2.
//
// For a tournament R on m vertices with adjacency matrix A (A_ab = 1 iff a -> b) and a new vertex v with out-set X:
//   c_3(R + v_X) - c_3(R) = x^T A (1 - x),   c_5(R + v_X) - c_5(R) = x^T A^3 (1 - x)        (x = indicator of X),
// because a cycle of length l through v is a path on l - 1 vertices of R from X to its complement, and a walk
// a -> b -> c -> d with a != d in a tournament is a path.  Modulo 2, c_3 = tr A^3 and c_5 = tr A^5 (a closed walk of
// length 3 or 5 in a tournament is a cycle).  So for K = {3}, {5}, {3,5} the cube sum
//   G_K(R) = sum_X (-1)^(c_K(R + v_X))
// is (-1)^(c_K(R)) times the Gauss sum of the F_2-quadratic form x^T W (1 + x), W = A, A^3, A + A^3, and is computed
// here by symplectic elimination (no enumeration of the 2^m out-sets).
//
// Modes:
//   k35cube census m            < gentourng output of order m : every class, labelled weights m!/|Aut|;
//                                 checks: trace formula against direct enumeration of 5-cycles, elimination against
//                                 direct summation over the 2^m out-sets (m <= 9), K = {3} against Theorem T.
//   k35cube mc m N seed         : N uniformly random labelled tournaments of order m.
// Further checks on every base (Theorem P of the note):
//   (P1) c_3 = C(q_1, 2) + floor(m/4) and c_3 + c_5 = C(q_2, 2) + floor(m/4) (mod 2), q_1 = |s|, q_2 = |As|, s = A1 mod 2;
//   (P2) the closed formula for G_{3,5}(R) in terms of K_0 = ker A^T and K_1 = {u : A^T u = 1} over F_2 and the weights
//        of their elements mod 4 (elements of K_0 have weight 0 or 1, elements of K_1 weight 0 or 3 mod 4):
//          m = 0 (4): 2^(m/2) |K_0|;                         m = 1 (4): 2^((m-1)/2) (|K_0| + n13 - n10);
//          m = 2 (4): 2^(m/2) (n00 + n13);                   m = 3 (4): 2^((m+1)/2) (n00 + n10).
// Output (one line per K): the exact signed labelled sum (census) or the Monte Carlo sums (mc), the numbers of bases
// with G_K < 0, = 0, > 0, and the histogram of the rank of the polar form.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef unsigned long long u64; typedef long long i64;
#define MAXM 62
static int m; static u64 row[MAXM], full;
static u64 sm_state;
static inline u64 splitmix(void){ u64 z=(sm_state+=0x9E3779B97F4A7C15ULL); z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL; z=(z^(z>>27))*0x94D049BB133111EBULL; return z^(z>>31); }
static inline int par(u64 x){ return __builtin_parityll(x); }
// rows of the product (as bit rows): P[a] = XOR over b in row[a] of M[b]
static void mulA(const u64 *M, u64 *P){ for(int a=0;a<m;a++){ u64 acc=0; for(u64 t=row[a];t;t&=t-1) acc^=M[__builtin_ctzll(t)]; P[a]=acc; } }
// Gauss sum of q(x) = sum_{a != d} W_ad x_a (1 + x_d) over F_2^m, divided by 2^m: returns 0 or +-2^(-rank/2).
static double gauss(const u64 *Wr, int *rank){
  u64 B[MAXM], Wt[MAXM]; u64 l=0;
  for(int a=0;a<m;a++) Wt[a]=0;
  for(int a=0;a<m;a++){ u64 w=Wr[a]&full&~(1ULL<<a); if(par(w)) l|=1ULL<<a; for(u64 t=w;t;t&=t-1) Wt[__builtin_ctzll(t)]|=1ULL<<a; }
  for(int a=0;a<m;a++) B[a]=((Wr[a]&full)^Wt[a])&~(1ULL<<a);
  u64 active=full; int c=0, r=0;
  for(;;){ int a=-1; for(u64 t=active;t;t&=t-1){ int i=__builtin_ctzll(t); if(B[i]&active){a=i;break;} }
    if(a<0) break; int d=__builtin_ctzll(B[a]&active);
    active&=~((1ULL<<a)|(1ULL<<d)); u64 Aa=B[a]&active, Bd=B[d]&active; int la=(l>>a)&1, ld=(l>>d)&1;
    for(u64 t=Aa;t;t&=t-1) B[__builtin_ctzll(t)]^=Bd;
    for(u64 t=Bd;t;t&=t-1) B[__builtin_ctzll(t)]^=Aa;
    l^=(Aa&Bd); if(la) l^=Bd; if(ld) l^=Aa; c^=(la&ld); r+=2; }
  *rank=r; if(l&active) return 0.0; double v=ldexp(1.0,-(r/2)); return c?-v:v; }
// direct summation (check)
static double direct(const u64 *Wr){ i64 tot=0; for(u64 X=0;X<=full;X++){ int q=0; for(u64 t=X;t;t&=t-1){ int a=__builtin_ctzll(t); q^=par(Wr[a]&~X&full&~(1ULL<<a)); } tot+= q?-1:1; } return (double)tot/(double)(full+1); }
// null space of the system [A^T | 1] (u, tau) = 0 over F_2: unknowns u_0..u_{m-1}, tau = bit m.  Counts the
// elements by tau and weight of u mod 4 into cnt[tau][w mod 4]; returns the dimension.
static int kernel_counts(u64 cnt[2][4]){
  u64 eq[MAXM]; for(int b=0;b<m;b++){ u64 r=0; for(int a=0;a<m;a++) if(row[a]>>b&1) r|=1ULL<<a; eq[b]=r|(1ULL<<m); }
  int piv[MAXM+1]; for(int j=0;j<=m;j++) piv[j]=-1; int rk=0;
  for(int j=0;j<=m && rk<m;j++){ int p=-1; for(int i=rk;i<m;i++) if(eq[i]>>j&1){p=i;break;} if(p<0) continue;
    u64 t=eq[p]; eq[p]=eq[rk]; eq[rk]=t; for(int i=0;i<m;i++) if(i!=rk && (eq[i]>>j&1)) eq[i]^=eq[rk]; piv[j]=rk; rk++; }
  u64 basis[MAXM+1]; int d=0;
  for(int f=0;f<=m;f++){ if(piv[f]>=0) continue; u64 v=1ULL<<f; for(int j=0;j<=m;j++) if(piv[j]>=0 && (eq[piv[j]]>>f&1)) v|=1ULL<<j; basis[d++]=v; }
  for(int t=0;t<2;t++)for(int w=0;w<4;w++) cnt[t][w]=0;
  if(d>24){ fprintf(stderr,"kernel too large\n"); exit(1); }
  u64 v=0; for(u64 g=0;g<(1ULL<<d);g++){ if(g){ int bit=__builtin_ctzll(g); v^=basis[bit]; } // Gray code
    int tau=(v>>m)&1; int w=__builtin_popcountll(v&full)&3; cnt[tau][w]++; }
  return d; }
static int perm[MAXM], sc[MAXM]; static u64 usedm, autcount;
static void autrec(int k){ if(k==m){autcount++;return;} for(int c=0;c<m;c++){ if(usedm>>c&1) continue; if(sc[c]!=sc[k]) continue; int ok=1;
    for(int j=0;j<k;j++){ int a=(row[k]>>j)&1, b=(row[c]>>perm[j])&1; if(a!=b){ok=0;break;} } if(!ok) continue; perm[k]=c; usedm|=1ULL<<c; autrec(k+1); usedm&=~(1ULL<<c); } }
static u64 c5direct(void){ u64 c5=0; for(int v=0;v<m;v++){ u64 hi=full&~((1ULL<<(v+1))-1);
    for(u64 a=row[v]&hi;a;a&=a-1){ int x=__builtin_ctzll(a); for(u64 b=row[x]&hi&~(1ULL<<x);b;b&=b-1){ int y=__builtin_ctzll(b);
      for(u64 c=row[y]&hi&~((1ULL<<x)|(1ULL<<y));c;c&=c-1){ int z=__builtin_ctzll(c);
        for(u64 d=row[z]&hi&~((1ULL<<x)|(1ULL<<y)|(1ULL<<z));d;d&=d-1){ int w=__builtin_ctzll(d); if(row[w]>>v&1) c5++; } } } } } return c5; }
int main(int argc,char**argv){
  int census=!strcmp(argv[1],"census"); m=atoi(argv[2]); full=(m==64)?~0ULL:((1ULL<<m)-1);
  u64 N=census?0:strtoull(argv[3],0,10); if(!census) sm_state=strtoull(argv[4],0,10)*0x2545F4914F6CDD1DULL+4242;
  static const char *nm[3]={"{3}","{5}","{3,5}"};
  double S1[3]={0,0,0}, S2[3]={0,0,0}; u64 neg[3]={0,0,0}, zer[3]={0,0,0}, pos[3]={0,0,0}; static u64 rk[3][MAXM+2];
  // census: exact integer sums of weight * G_K (G_K = b * 2^m is an integer)
  __int128 ex[3]={0,0,0}; u64 wtot=0, nt=0, chk5=0, chkd=0, chkT=0, chkP1=0, chkP2=0, chkW=0; u64 fact=1; for(int i=2;i<=m;i++) fact*=i;
  char line[512]; u64 it=0;
  for(;;){
    if(census){ if(!fgets(line,sizeof line,stdin)) break; if(line[0]!='0'&&line[0]!='1') continue; int k=0;
      for(int i=0;i<m;i++) row[i]=0;
      for(int i=0;i<m;i++)for(int j=i+1;j<m;j++){ if(line[k++]=='1') row[i]|=1ULL<<j; else row[j]|=1ULL<<i; } }
    else { if(it++>=N) break; for(int i=0;i<m;i++) row[i]=0; u64 bits=0; int nb=0;
      for(int i=0;i<m;i++)for(int j=i+1;j<m;j++){ if(nb==0){bits=splitmix();nb=64;} if(bits&1) row[i]|=1ULL<<j; else row[j]|=1ULL<<i; bits>>=1; nb--; } }
    nt++;
    u64 A2[MAXM],A3[MAXM],A4[MAXM],A5[MAXM],W[MAXM];
    mulA(row,A2); mulA(A2,A3); mulA(A3,A4); mulA(A4,A5);
    int t3=0,t5=0; for(int a=0;a<m;a++){ t3^=(A3[a]>>a)&1; t5^=(A5[a]>>a)&1; }
    { // (P1)
      u64 svec=0; for(int a=0;a<m;a++) if(__builtin_popcountll(row[a])&1) svec|=1ULL<<a;
      u64 As=0; for(int a=0;a<m;a++) if(par(row[a]&svec)) As|=1ULL<<a;
      int q1=__builtin_popcountll(svec), q2=__builtin_popcountll(As), kap=(m/4)&1;
      if((((q1*(q1-1)/2)+kap)&1)!=t3) chkP1++;
      if((((q2*(q2-1)/2)+kap)&1)!=(t3^t5)) chkP1++; }
    u64 w=1;
    if(census){ if((c5direct()&1)!=(u64)t5) chk5++;
      i64 c3=(i64)m*(m-1)*(m-2)/6; for(int i=0;i<m;i++){ sc[i]=__builtin_popcountll(row[i]); c3-=(i64)sc[i]*(sc[i]-1)/2; } if((c3&1)!=t3) chk5++;
      autcount=0; usedm=0; autrec(0); w=fact/autcount; wtot+=w; }
    for(int K=0;K<3;K++){ int eps = (K==0)?t3:(K==1)?t5:(t3^t5);
      for(int a=0;a<m;a++) W[a]= (K==0)?row[a]:(K==1)?A3[a]:(row[a]^A3[a]);
      int r; double b=gauss(W,&r); if(eps) b=-b;
      if(census && m<=9){ double d=direct(W); if(eps) d=-d; if(fabs(d-b)>1e-12) chkd++; }
      if(K==0){ double T=ldexp(1.0,(m+1)/2-m); if(fabs(b-T)>1e-15) chkT++; }
      if(K==2){ u64 kc[2][4]; kernel_counts(kc);
        if(kc[0][2]||kc[0][3]||kc[1][1]||kc[1][2]) chkW++;
        double F; double k0=(double)(kc[0][0]+kc[0][1]);
        switch(m&3){ case 0: F=ldexp(k0,m/2); break;
                     case 1: F=ldexp(k0+(double)kc[1][3]-(double)kc[1][0],(m-1)/2); break;
                     case 2: F=ldexp((double)(kc[0][0]+kc[1][3]),m/2); break;
                     default: F=ldexp((double)(kc[0][0]+kc[1][0]),(m+1)/2); }
        if(fabs(ldexp(b,m)-F)>1e-6*fabs(F)+1e-9) chkP2++; }
      rk[K][r]++; if(b<0) neg[K]+=1; else if(b==0) zer[K]+=1; else pos[K]+=1;
      if(census){ // G = b 2^m, an integer
        i64 Gi=(i64)llround(ldexp(b,m)); ex[K]+=(__int128)Gi*(__int128)w; }
      else { S1[K]+=b; S2[K]+=b*b; } }
  }
  for(int K=0;K<3;K++){
    if(census){ __int128 e=ex[K]; int sg=e<0; if(sg) e=-e; char buf[64]; int p=63; buf[p]=0; if(e==0) buf[--p]='0'; while(e){ buf[--p]='0'+(int)(e%10); e/=10; }
      printf("census m %d K %s classes %llu weight %llu signedsum %s%s neg %llu zero %llu pos %llu rank", m, nm[K], nt, wtot, sg?"-":"", buf+p, neg[K], zer[K], pos[K]); }
    else printf("mc m %d K %s N %llu S1 %.17g S2 %.17g neg %llu zero %llu pos %llu rank", m, nm[K], nt, S1[K], S2[K], neg[K], zer[K], pos[K]);
    for(int r=0;r<=m;r+=2) if(rk[K][r]) printf(" r%d %llu", r, rk[K][r]);
    printf("\n"); }
  fprintf(stderr,"m %d: bases %llu; check failures: trace formulas %llu, direct summation %llu, Theorem T %llu, (P1) score formulas %llu, (P2) closed formula %llu, kernel weights %llu\n", m, nt, chk5, chkd, chkT, chkP1, chkP2, chkW);
  return 0; }
