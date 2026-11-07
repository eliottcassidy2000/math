/* resdp.c - kind-pasteur-2026-11-07-S15. First-principles local factors of the two-root-change word multiplicity.
 * The moves O(D,E,U,d) = (D,0,2U,d), Z = (D-E+U,d,2U-d,d), J = (D-E+U,D-E,D+U,D) are integer-linear, so the residues
 * (D,E,U,d) mod q evolve by a finite automaton. Count every prefix w with at most KMAX letters J and exactly n letters
 * O/Z (J letters anywhere, also at the end), with weight 2 if #J(w) < KMAX (it is counted as wZ and wJ in M_K) and 1 if
 * #J(w) = KMAX (only wZ), and record the class of its value F(t.w) = D - E + U (mod q). The local factor at stage n is
 * L_n(r) = q * (weight in class r) / (total weight), averaged over the odd N (the q classes are equinumerous among
 * odd N because q is odd). No census data are used.
 * usage: resdp q nmax KMAX   ; prints L_n by class aggregate for n = 5, 10, ..., nmax
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
static int q, KMAX; static ll Q4;
static inline ll enc(int D,int E,int U,int d){ return (((ll)D*q + E)*q + U)*q + d; }
static inline void dec(ll s, int*D,int*E,int*U,int*d){ *d = s%q; s/=q; *U = s%q; s/=q; *E = s%q; s/=q; *D = (int)s; }
static inline int md(ll x){ x%=q; if(x<0) x+=q; return (int)x; }
static int vq(int r, int p){ if (r==0) return 99; int k=0; while (r%p==0){ r/=p; k++; } return k; }
int main(int argc, char**argv){
  if (argc<4){ fprintf(stderr,"usage: resdp q nmax KMAX\n"); return 1; }
  q = atoi(argv[1]); int nmax = atoi(argv[2]); KMAX = atoi(argv[3]);
  Q4 = (ll)q*q*q*q;
  double **cur = malloc((KMAX+1)*sizeof(double*)), **nxt = malloc((KMAX+1)*sizeof(double*));
  for (int j=0;j<=KMAX;j++){ cur[j]=calloc(Q4,sizeof(double)); nxt[j]=calloc(Q4,sizeof(double)); }
  cur[0][enc(1%q,0,2%q,1%q)] = 1.0;           /* t = (1,0,2,1) */
  int p = q; for (int f=2; f<=q; f++) if (q%f==0){ p=f; break; }   /* smallest prime factor of q */
  for (int n=0; n<=nmax; n++){
    /* close under J (at most KMAX in total), stage n */
    for (int j=0;j<KMAX;j++) for (ll s=0;s<Q4;s++){ double w=cur[j][s]; if(!w) continue; int D,E,U,d; dec(s,&D,&E,&U,&d);
      cur[j+1][enc(md(D-E+U), md(D-E), md(D+U), D)] += w; }
    if (n%5==0 || n==nmax){
      double *cls = calloc(q,sizeof(double)); double tot=0;
      for (int j=0;j<=KMAX;j++){ double wt = (j<KMAX)? 2.0 : 1.0;
        for (ll s=0;s<Q4;s++){ double w=cur[j][s]; if(!w) continue; int D,E,U,d; dec(s,&D,&E,&U,&d); cls[md(D-E+U)] += wt*w; tot += wt*w; } }
      /* aggregate by p-adic valuation of the class (r = 0 means v >= log_p q) */
      int kq=0; { int t=q; while(t%p==0){ t/=p; kq++; } }
      printf("q=%d n=%d", q, n);
      double sv[32]={0}, nv[32]={0};
      for (int r=0;r<q;r++){ int v = vq(r,p); if (v>kq) v=kq; sv[v]+=cls[r]; nv[v]+=1; }
      printf("  L by v_%d(N)=0..%d:", p, kq); for (int v=0; v<=kq; v++) printf(" %.4f", nv[v]? (sv[v]/nv[v])*q/tot : 0.0);
      if (q<=31){ printf("  | by class:"); for (int r=0;r<q;r++) printf(" %.4f", cls[r]*q/tot); }
      printf("\n"); free(cls);
    }
    if (n==nmax) break;
    /* one letter O or Z at every level */
    for (int j=0;j<=KMAX;j++) memset(nxt[j],0,Q4*sizeof(double));
    for (int j=0;j<=KMAX;j++) for (ll s=0;s<Q4;s++){ double w=cur[j][s]; if(!w) continue; int D,E,U,d; dec(s,&D,&E,&U,&d);
      nxt[j][enc(D, 0, md(2*U), d)] += w;                       /* O */
      nxt[j][enc(md(D-E+U), d, md(2*U-d), d)] += w; }           /* Z */
    /* renormalise to avoid overflow */
    double tot=0; for (int j=0;j<=KMAX;j++) for (ll s=0;s<Q4;s++) tot+=nxt[j][s];
    for (int j=0;j<=KMAX;j++){ for (ll s=0;s<Q4;s++) nxt[j][s]/=tot; double*t=cur[j]; cur[j]=nxt[j]; nxt[j]=t; }
  }
  return 0;
}
