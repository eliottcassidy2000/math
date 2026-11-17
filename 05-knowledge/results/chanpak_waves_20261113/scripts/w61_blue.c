/* Wave 61. */
/* Blue lines of the tiling-line metagraph (THM-4731/4733): tilings t of the tiles (i,j), j >= i+2, along the path 0->1->...->n-1
   (bit 0: i->j, bit 1: j->i); kappa = complement of all tiles; sigma = mirror (i,j) -> (n-1-j, n-1-i).  A line {t, kappa t} is
   blue iff sigma fixes it; for n >= 4 this means sigma t = t.  BFS in the blue graph on classes from the transitive class. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static int n, m, ti[64], tj[64], mir[64];
static void build(u64 t, unsigned*out){ for(int v=0;v<n;v++) out[v]=0; for(int v=0;v+1<n;v++) out[v]|=1u<<(v+1);
  for(int k=0;k<m;k++){ int i=ti[k], j=tj[k]; if(t>>k&1) out[j]|=1u<<i; else out[i]|=1u<<j; } }
static u64 Hcount(const unsigned*o){ static u64 dp[1<<10][10]; for(int S=0;S<(1<<n);S++) for(int v=0;v<n;v++) dp[S][v]=0; for(int v=0;v<n;v++) dp[1<<v][v]=1;
  for(int S=1;S<(1<<n);S++) for(int v=0;v<n;v++) if(dp[S][v]){ unsigned t=o[v]&~S; while(t){ int w=__builtin_ctz(t); t&=t-1; dp[S|1u<<w][w]+=dp[S][v]; } }
  u64 h=0; for(int v=0;v<n;v++) h+=dp[(1<<n)-1][v]; return h; }
/* canonical code: WL colour refinement then brute force within colour classes */
static int col[10], perm[10], used[10]; static u64 best; static const unsigned*O;
static u64 code(void){ u64 c=0; int b=0; for(int i=0;i<n;i++) for(int j=i+1;j<n;j++,b++) if(O[perm[i]]>>perm[j]&1) c|=1ULL<<b; return c; }
static int order_[10];
static void rec(int pos){ if(pos==n){ u64 c=code(); if(c<best) best=c; return; }
  int target=col[order_[pos]]; for(int v=0;v<n;v++){ if(used[v]||col[v]!=target) continue; used[v]=1; perm[pos]=v; rec(pos+1); used[v]=0; } }
static u64 canon(const unsigned*o){ O=o; for(int v=0;v<n;v++) col[v]=__builtin_popcount(o[v]);
  for(int it=0;it<n;it++){ long long sig[10]; for(int v=0;v<n;v++){ long long s=col[v]*1000003LL; for(int w=0;w<n;w++) if(w!=v){ if(o[v]>>w&1) s+= (long long)(col[w]+1)*7919; else s+=(long long)(col[w]+1)*104729; } sig[v]=s; }
    /* relabel signatures to small ints preserving order */
    int nc[10]; for(int v=0;v<n;v++){ int r=0; for(int w=0;w<n;w++) if(sig[w]<sig[v]) r++; nc[v]=r; } int changed=0; for(int v=0;v<n;v++) if(nc[v]!=col[v]) changed=1; memcpy(col,nc,sizeof nc); if(!changed) break; }
  /* order positions by colour */ int k=0; for(int c=0;c<n;c++) for(int v=0;v<n;v++) if(col[v]==c) order_[k++]=v;
  best=~0ULL; memset(used,0,sizeof used); rec(0); return best; }
typedef struct { u64 key; int id; } KV;
int main(int argc,char**argv){ int N1=atoi(argv[1]);
  for(n=3;n<=N1;n++){ m=0; for(int i=0;i<n;i++) for(int j=i+2;j<n;j++){ ti[m]=i; tj[m]=j; m++; }
    for(int k=0;k<m;k++){ int i2=n-1-tj[k], j2=n-1-ti[k]; for(int q=0;q<m;q++) if(ti[q]==i2&&tj[q]==j2) mir[k]=q; }
    /* enumerate sigma-fixed tilings: free bits = one per sigma-orbit of tiles */
    int orb[64], no=0, rep[64]; for(int k=0;k<m;k++){ rep[k]=-1; } for(int k=0;k<m;k++){ if(rep[k]>=0) continue; rep[k]=no; rep[mir[k]]=no; orb[no++]=k; }
    long nfix=1L<<no; /* classes */ int cap=200000; u64*keys=malloc(sizeof(u64)*cap); u64*Hval=malloc(sizeof(u64)*cap); int ncls=0;
    int *ea=malloc(sizeof(int)*nfix), *eb=malloc(sizeof(int)*nfix); long ne=0; int tclass=-1;
    for(long f=0;f<nfix;f++){ u64 t=0; for(int k=0;k<m;k++) if(f>>rep[k]&1) t|=1ULL<<k; u64 kt=t^((m==64)?~0ULL:((1ULL<<m)-1));
      int ids[2]; for(int s=0;s<2;s++){ unsigned out[10]; build(s?kt:t,out); u64 c=canon(out); int id=-1; for(int q=0;q<ncls;q++) if(keys[q]==c){ id=q; break; }
        if(id<0){ id=ncls; keys[ncls]=c; Hval[ncls]=Hcount(out); ncls++; } ids[s]=id; }
      if(t==0) tclass=ids[0];
      ea[ne]=ids[0]; eb[ne]=ids[1]; ne++; }
    /* BFS */ int *dist=malloc(sizeof(int)*ncls); for(int q=0;q<ncls;q++) dist[q]=-1; dist[tclass]=0; int changed=1, maxd=0;
    while(changed){ changed=0; for(long e=0;e<ne;e++){ int a=ea[e], b=eb[e]; if(dist[a]>=0&&dist[b]<0){ dist[b]=dist[a]+1; changed=1; if(dist[b]>maxd) maxd=dist[b]; } if(dist[b]>=0&&dist[a]<0){ dist[a]=dist[b]+1; changed=1; if(dist[a]>maxd) maxd=dist[a]; } } }
    int reach=0; for(int q=0;q<ncls;q++) if(dist[q]>=0) reach++;
    printf("n=%d: sigma-fixed tilings %ld; self-converse classes on blue lines %d, reached from the transitive class %d; layers:",n,nfix,ncls,reach);
    for(int d=0;d<=maxd;d++){ int c=0; for(int q=0;q<ncls;q++) if(dist[q]==d) c++; printf(" %d",c); } printf("\n");
    printf("   H at distance 1:"); for(int q=0;q<ncls;q++) if(dist[q]==1) printf(" %llu",Hval[q]); printf("   | distance 2:"); int c2=0; for(int q=0;q<ncls&&c2<16;q++) if(dist[q]==2){ printf(" %llu",Hval[q]); c2++; } printf("\n");
    /* degree of the transitive class (number of blue lines at it, counting multiplicity) */ long deg=0; for(long e=0;e<ne;e++) if(ea[e]==tclass||eb[e]==tclass) deg++; printf("   blue lines at the transitive class (tiling pairs): %ld\n",deg/1);
    free(keys); free(Hval); free(ea); free(eb); free(dist); }
  return 0; }
