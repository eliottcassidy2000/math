/* BFS of SEAM states (H,S,R) under a list of 3x3 gadget matrices; states with H<=X, R<=X expanded. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
typedef struct {ll h,s,r;} st;
static st *tab; static size_t cap; static size_t cnt;
static size_t hsh(ll h,ll s,ll r){ unsigned long long x=(unsigned long long)h*0x9E3779B97F4A7C15ULL ^ (unsigned long long)s*0xC2B2AE3D27D4EB4FULL ^ (unsigned long long)r*0x165667B19E3779F9ULL; x^=x>>29; return (size_t)(x&(cap-1)); }
static int ins(ll h,ll s,ll r){ size_t i=hsh(h,s,r); while(tab[i].h){ if(tab[i].h==h&&tab[i].s==s&&tab[i].r==r) return 0; i=(i+1)&(cap-1);} tab[i].h=h;tab[i].s=s;tab[i].r=r; cnt++; return 1; }
int main(int argc,char**argv){
  ll X=atoll(argv[1]); const char*mf=argv[2]; int onlyP = argc>3 ? atoi(argv[3]) : 0;
  FILE*f=fopen(mf,"r"); ll (*M)[9]=malloc(sizeof(ll)*9*40000); int nm=0;
  while(nm<40000 && fscanf(f,"%lld %lld %lld %lld %lld %lld %lld %lld %lld",&M[nm][0],&M[nm][1],&M[nm][2],&M[nm][3],&M[nm][4],&M[nm][5],&M[nm][6],&M[nm][7],&M[nm][8])==9) nm++;
  fclose(f);
  cap=1ULL<<26; tab=calloc(cap,sizeof(st));
  st *q=malloc(sizeof(st)*(cap/2)); size_t qh=0,qt=0;
  /* seeds: realisable triples of order <= 5 */
  ll seeds[][3]={{1,1,1},{1,0,2},{3,1,3},{1,0,4},{1,0,8},{1,0,16},{1,1,1}};
  char *hit=calloc(X+1,1);
  for(int i=0;i<6;i++){ if(ins(seeds[i][0],seeds[i][1],seeds[i][2])){ q[qt++]=(st){seeds[i][0],seeds[i][1],seeds[i][2]}; hit[seeds[i][0]]=1; } }
  while(qh<qt){ st a=q[qh++]; if(a.r> X) continue;
    for(int m=0;m<nm;m++){ ll *A=M[m];
      ll h=A[0]*a.h+A[1]*a.s+A[2]*a.r; if(h>X||h<=0) continue;
      ll s=A[3]*a.h+A[4]*a.s+A[5]*a.r; ll r=A[6]*a.h+A[7]*a.s+A[8]*a.r;
      if(onlyP && s!=0) continue;
      hit[h]=1;
      if(r<=X && ins(h,s,r)){ if(qt>=cap/2){fprintf(stderr,"queue full\n");return 1;} q[qt++]=(st){h,s,r}; }
    }
  }
  ll miss=0; printf("X=%lld mats=%d states=%zu\nmissing odd:",X,nm,cnt);
  for(ll n=1;n<=X;n+=2) if(!hit[n]){ miss++; if(miss<=60) printf(" %lld",n);} printf("\ncount missing=%lld\n",miss);
  /* dump P-plane states for fill analysis */
  FILE*o=fopen("pstates.txt","w"); for(size_t i=0;i<cap;i++) if(tab[i].h && tab[i].s==0) fprintf(o,"%lld %lld\n",tab[i].h,tab[i].r); fclose(o);
  return 0;
}
