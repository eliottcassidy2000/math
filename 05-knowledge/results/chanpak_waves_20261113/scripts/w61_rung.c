/* Wave 61. Usage: w61_rung NMAX X STORE MAXQ   (reads w57_reps.txt) */
/* k-smooth tournament profiles by recursion (C port of w61_rung.py).  Profile p[d] (d = 1..n) = ordered covers of V by d nonempty
   paths; chi(t) = sum_d C(t,d) p[d].  Join: chi multiplies.  Substitution Q[S_1..S_q]: chi(t) = sum_r N_t(r) prod p_i(r_i) with
   N_t(r) = coefficient of x^r in (1 + W_Q(x))^t, W_Q = walk generating function of Q (THM-3121).  Store blocks with H <= STORE,
   record counts of strong (substituted) tournaments and joins with H <= X.  Quotients: prime tournaments of order 3, 5, 6. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll; typedef __int128 i128;
#define MAXN 14
typedef struct { int n; ll p[MAXN+1]; } Prof;
static Prof *P[MAXN+1]; static int NP[MAXN+1], CP[MAXN+1]; static Prof *SP[MAXN+1]; static int NS[MAXN+1], CS[MAXN+1];
static ll binom[MAXN+2][MAXN+2];
static int NMAX; static ll X, STORE; static unsigned char *spec;
/* hash sets per order */
#define HS (1<<20)
static unsigned long long *hk[MAXN+1];
static unsigned long long phash(const Prof*a){ unsigned long long h=1469598103934665603ULL^a->n; for(int d=1;d<=a->n;d++){ h^=(unsigned long long)a->p[d]; h*=1099511628211ULL; } return h?h:1; }
static int addprof(Prof**arr,int*cnt,int*cap,unsigned long long*tab,const Prof*a){ unsigned long long h=phash(a); unsigned long long i=h&(HS-1);
  while(tab[i]){ if(tab[i]==h) return 0; i=(i+1)&(HS-1); } tab[i]=h; if(*cnt==*cap){ *cap=*cap?*cap*2:64; *arr=realloc(*arr,sizeof(Prof)**cap); } (*arr)[(*cnt)++]=*a; return 1; }
static void chi_of(const Prof*a,int T,i128*ch){ for(int t=0;t<=T;t++){ i128 s=0; for(int d=1;d<=a->n&&d<=t;d++) s+=(i128)binom[t][d]*a->p[d]; ch[t]=s; } }
static void p_of_chi(const i128*ch,int n,Prof*out){ out->n=n; for(int d=1;d<=n;d++){ i128 s=0; for(int t=0;t<=d;t++){ i128 term=(i128)binom[d][t]*ch[t]; if((d-t)&1) s-=term; else s+=term; } out->p[d]=(ll)s; } }
/* quotients */
static int nq, qord[32], qout[32][8]; static ll qH[32];
/* walk kernel for (Q, sizes): enumerate visit vectors in mixed radix */
static int rad[8], q_cur; static int nst;
static ll *Wpoly; /* W[idx] counts nonempty walks with visit vector idx */
static void walks(int qi,const int*sz){ int q=qord[qi]; nst=1; for(int i=0;i<q;i++){ rad[i]=sz[i]+1; nst*=rad[i]; }
  Wpoly=calloc(nst,sizeof(ll)); /* DP over (current vertex, visit vector) */
  ll *cur=calloc((size_t)nst*q,sizeof(ll)), *nxt=calloc((size_t)nst*q,sizeof(ll));
  int mul[8]; mul[0]=1; for(int i=1;i<q;i++) mul[i]=mul[i-1]*rad[i-1];
  for(int v=0;v<q;v++) if(sz[v]>=1) cur[(size_t)mul[v]*q+v]+=1;
  int any=1; while(any){ any=0; memset(nxt,0,sizeof(ll)*(size_t)nst*q);
    for(int idx=0;idx<nst;idx++) for(int v=0;v<q;v++){ ll c=cur[(size_t)idx*q+v]; if(!c) continue; Wpoly[idx]+=c;
      int vis[8]; int r=idx; for(int i=0;i<q;i++){ vis[i]=r%rad[i]; r/=rad[i]; }
      for(int w=0;w<q;w++){ if(w==v||!(qout[qi][v]>>w&1)) continue; if(vis[w]+1>sz[w]) continue; nxt[(size_t)(idx+mul[w])*q+w]+=c; any=1; } }
    ll*tmp=cur; cur=nxt; nxt=tmp; }
  free(cur); free(nxt); }
static void compose_q(int qi,const int*sz,int n){ int q=qord[qi]; walks(qi,sz); int mul[8]; mul[0]=1; for(int i=1;i<q;i++) mul[i]=mul[i-1]*rad[i-1];
  /* N_t for t = 0..n */
  ll **N=malloc(sizeof(ll*)*(n+1)); N[0]=calloc(nst,sizeof(ll)); N[0][0]=1;
  for(int t=1;t<=n;t++){ N[t]=calloc(nst,sizeof(ll)); for(int a=0;a<nst;a++){ ll c=N[t-1][a]; if(!c) continue; N[t][a]+=c; int va[8]; int r=a; for(int i=0;i<q;i++){ va[i]=r%rad[i]; r/=rad[i]; }
        for(int b=1;b<nst;b++){ ll w=Wpoly[b]; if(!w) continue; int vb[8]; int r2=b, ok=1; for(int i=0;i<q;i++){ vb[i]=r2%rad[i]; r2/=rad[i]; if(va[i]+vb[i]>sz[i]){ ok=0; break; } } if(!ok) continue; N[t][a+b]+=c*w; } } }
  /* list of r vectors with all r_i >= 1 */
  int *rl=malloc(sizeof(int)*nst), nr=0; int (*rv)[8]=malloc(sizeof(int[8])*nst);
  for(int a=0;a<nst;a++){ int r=a, ok=1; for(int i=0;i<q;i++){ rv[nr][i]=r%rad[i]; r/=rad[i]; if(rv[nr][i]==0) ok=0; } if(ok) rl[nr++]=a; }
  /* iterate block tuples */
  int idx[8]={0}; Prof *lists[8]; int lens[8]; for(int i=0;i<q;i++){ lists[i]=P[sz[i]]; lens[i]=NP[sz[i]]; if(!lens[i]){ goto done; } }
  while(1){ i128 lb=qH[qi]; for(int i=0;i<q;i++) lb*=lists[i][idx[i]].p[1];
    if(lb<=X){ i128 ch[MAXN+2]; for(int t=0;t<=n;t++){ i128 s=0; for(int k=0;k<nr;k++){ int a=rl[k]; ll c=N[t][a]; if(!c) continue; i128 pr=c; for(int i=0;i<q;i++) pr*=lists[i][idx[i]].p[rv[k][i]]; s+=pr; } ch[t]=s; }
      Prof out; p_of_chi(ch,n,&out); if(out.p[1]<=X){ spec[out.p[1]]=1; }
      if(out.p[1]<=STORE){ addprof(&SP[n],&NS[n],&CS[n],hk[0],&out); } }
    int i=0; while(i<q){ idx[i]++; if(idx[i]<lens[i]) break; idx[i]=0; i++; } if(i==q) break; }
  done: for(int t=0;t<=n;t++) free(N[t]); free(N); free(rl); free(rv); free(Wpoly); }
static void comps(int n,int q,int pos,int*sz,int qi){ if(pos==q-1){ sz[pos]=n; if(n>=1) compose_q(qi,sz,0+(sz[0]*0)+ (int)0 + n + 0 - n + ({int s=0; for(int i=0;i<q;i++) s+=sz[i]; s;}) ); return; }
  for(int a=1;a<=n-(q-1-pos);a++){ sz[pos]=a; comps(n-a,q,pos+1,sz,qi); } }
int main(int argc,char**argv){ NMAX=atoi(argv[1]); X=atoll(argv[2]); STORE=atoll(argv[3]); int maxq=atoi(argv[4]);
  for(int i=0;i<=MAXN+1;i++){ binom[i][0]=1; for(int j=1;j<=i;j++) binom[i][j]=binom[i-1][j-1]+(j<=i-1?binom[i-1][j]:0); }
  spec=calloc(X+1,1);
  /* primes from w57_reps.txt */ FILE*fi=fopen("w57_reps.txt","r"); char line[256];
  while(fgets(line,sizeof line,fi)){ char*p=line; int n=strtol(p,&p,10); int o[8]; for(int i=0;i<n;i++) o[i]=strtol(p,&p,10); if(n!=3&&n!=5&&n!=6) continue; if(n>maxq) continue;
    int prime=1; for(int M=1;M<(1<<n)-1&&prime;M++){ if(__builtin_popcount(M)<2) continue; int ok=1; for(int x=0;x<n&&ok;x++) if(!(M>>x&1)){ int a=-1; for(int y=0;y<n;y++) if(M>>y&1){ int b=o[x]>>y&1; if(a<0) a=b; else if(a!=b){ ok=0; break; } } } if(ok) prime=0; }
    if(!prime) continue; qord[nq]=n; memcpy(qout[nq],o,sizeof(int)*n);
    ll h=0; int perm[8]; for(int i=0;i<n;i++) perm[i]=i; /* H by brute force */ { int c[8]={0}; int k=0; while(1){ int good=1; for(int i=0;i+1<n;i++) if(!(o[perm[i]]>>perm[i+1]&1)){ good=0; break; } h+=good; if(k>=n) break; while(k<n){ if(c[k]<k){ if(k%2==0){int t=perm[0];perm[0]=perm[k];perm[k]=t;} else {int t=perm[c[k]];perm[c[k]]=perm[k];perm[k]=t;} c[k]++; k=0; break; } else { c[k]=0; k++; } } if(k>=n) break; } }
    qH[nq]=h; nq++; }
  fclose(fi); printf("prime quotients loaded: %d\n",nq);
  for(int i=0;i<=MAXN;i++) hk[i]=calloc(HS,sizeof(unsigned long long));
  Prof one; one.n=1; one.p[1]=1; addprof(&P[1],&NP[1],&CP[1],hk[1],&one); spec[1]=1; NS[1]=0; addprof(&SP[1],&NS[1],&CS[1],hk[0],&one);
  for(int n=2;n<=NMAX;n++){ memset(hk[0],0,sizeof(unsigned long long)*HS);
    for(int qi=0;qi<nq;qi++){ int q=qord[qi]; if(q>n) continue; int sz[8]; comps(n,q,0,sz,qi); }
    /* stored strong at order n are SP[n]; all stored = strong plus joins (first component strong) */
    for(int k=0;k<NS[n];k++) addprof(&P[n],&NP[n],&CP[n],hk[n],&SP[n][k]);
    for(int a=1;a<n;a++) for(int i=0;i<NS[a];i++) for(int j=0;j<NP[n-a];j++){ ll h=SP[a][i].p[1]*P[n-a][j].p[1]; if(h<=X) spec[h]=1; if(h>STORE) continue;
        i128 ca[MAXN+2], cb[MAXN+2], cc[MAXN+2]; chi_of(&SP[a][i],n,ca); chi_of(&P[n-a][j],n,cb); for(int t=0;t<=n;t++) cc[t]=ca[t]*cb[t]; Prof out; p_of_chi(cc,n,&out); addprof(&P[n],&NP[n],&CP[n],hk[n],&out); }
    long rec=0; for(ll h=1;h<=X;h+=2) if(spec[h]) rec++;
    printf("order %d: stored %d (strong %d); recorded odd counts <= %lld: %ld\n",n,NP[n],NS[n],X,rec); fflush(stdout); }
  printf("odd counts <= %lld NOT realized:",X); for(ll h=1;h<=X;h+=2) if(!spec[h]) printf(" %lld",h); printf("\n"); return 0; }
