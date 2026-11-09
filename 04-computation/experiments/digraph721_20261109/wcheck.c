/* wcheck.c - independent brute-force check of explicit witness graphs (THM-5055 base data).
   Input lines: "tag set n k | u-v u-v ...".  For each graph: d(G) by subset enumeration must equal k;
   von Bromssen parity via full stabiliser-chain transversals (generators of Aut(G)) must be EVEN;
   the flags required by the set must hold:
     W: no universal vertex; P: a pendant vertex; S: a pendant whose neighbour is non-universal;
     S2: >= 2 such pendants; S3: a vertex of degree 2; G: S2 and S3. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define MAXN 32
typedef uint32_t row;
static int n; static row adj[MAXN];
static uint32_t *dom;
static uint64_t count_dom(void){
  uint32_t full = (n==32)?0xffffffffu:((1u<<n)-1); uint64_t c=0; uint64_t N=1ull<<n;
  row cl[MAXN]; for(int v=0;v<n;v++) cl[v]=adj[v]|(1u<<v);
  dom[0]=0;
  for(uint64_t S=1;S<N;S++){ int b=__builtin_ctzll(S); dom[S]=dom[S&(S-1)]|cl[b]; if(dom[S]==full) c++; }
  return c;
}
/* colour refinement (exact partition refinement, not hashed) */
static int col[MAXN];
static void refine(void){
  int c[MAXN], nc[MAXN]; for(int v=0;v<n;v++) c[v]=__builtin_popcount(adj[v]);
  for(int round=0;round<n;round++){
    /* signature = (c[v], sorted multiset of neighbour colours) -> rank */
    static int sig[MAXN][MAXN+1];
    for(int v=0;v<n;v++){ int k=0; sig[v][k++]=c[v]; int tmp[MAXN],t=0; row r=adj[v]; while(r){int u=__builtin_ctz(r); r&=r-1; tmp[t++]=c[u];}
      for(int i=1;i<t;i++){int x=tmp[i],j=i-1; while(j>=0&&tmp[j]>x){tmp[j+1]=tmp[j];j--;} tmp[j+1]=x;}
      for(int i=0;i<t;i++) sig[v][k++]=tmp[i]; for(;k<=MAXN;k++) sig[v][k]=-1; }
    int ord[MAXN]; for(int v=0;v<n;v++) ord[v]=v;
    for(int i=1;i<n;i++){int x=ord[i],j=i-1; while(j>=0 && memcmp(sig[ord[j]],sig[x],sizeof(sig[0]))>0){ord[j+1]=ord[j];j--;} ord[j+1]=x;}
    int cur=0; nc[ord[0]]=0; for(int i=1;i<n;i++){ if(memcmp(sig[ord[i]],sig[ord[i-1]],sizeof(sig[0]))!=0) cur++; nc[ord[i]]=cur; }
    int same=1; /* stable if number of classes unchanged */
    int k1=0,k2=0; { int seen[MAXN*MAXN+2]; memset(seen,0,sizeof(seen)); }
    int distinct_old=0, distinct_new=cur+1; { int s[MAXN],t=0; for(int v=0;v<n;v++){int f=0; for(int i=0;i<t;i++) if(s[i]==c[v]){f=1;break;} if(!f) s[t++]=c[v];} distinct_old=t; }
    (void)k1;(void)k2;(void)same;
    for(int v=0;v<n;v++) c[v]=nc[v];
    if(distinct_new==distinct_old && round>0) break;
  }
  for(int v=0;v<n;v++) col[v]=c[v];
}
static int amap[MAXN]; static uint32_t aused;
static int aut_rec(int v){
  if(v==n) return 1;
  if(amap[v]>=0){ int c=amap[v]; for(int j=0;j<v;j++) if(((adj[v]>>j)&1)!=((adj[c]>>amap[j])&1)) return 0; return aut_rec(v+1); }
  for(int c=0;c<n;c++){ if((aused>>c)&1) continue; if(col[c]!=col[v]) continue;
    int ok=1; for(int j=0;j<v;j++) if(((adj[v]>>j)&1)!=((adj[c]>>amap[j])&1)){ok=0;break;}
    if(!ok) continue;
    amap[v]=c; aused|=1u<<c; if(aut_rec(v+1)) return 1; amap[v]=-1; aused&=~(1u<<c); }
  return 0;
}
static int sign_of(const int*p){ int par=0; for(int u=0;u<n;u++){ row r=adj[u]&~((u==31)?0xffffffffu:((1u<<(u+1))-1)); while(r){int v=__builtin_ctz(r); r&=r-1; if(p[u]>p[v]) par^=1;} } return par; }
static int is_odd(void){
  refine();
  for(int i=0;i<n;i++) for(int c=i+1;c<n;c++){ if(col[c]!=col[i]) continue;
    for(int v=0;v<n;v++) amap[v]=-1; aused=0;
    for(int j=0;j<i;j++){amap[j]=j; aused|=1u<<j;}
    amap[i]=c; aused|=1u<<c;
    if(aut_rec(0) && sign_of(amap)) return 1; }
  return 0;
}
int main(void){
  dom=malloc(sizeof(uint32_t)<<22);
  char line[1<<16]; long ok=0,bad=0;
  while(fgets(line,sizeof line,stdin)){
    char tag[64],set[8]; unsigned long long k; int pos;
    if(sscanf(line,"%63s %7s %d %llu |%n",tag,set,&n,&k,&pos)<4) continue;
    memset(adj,0,sizeof adj);
    char*p=line+pos; int u,v,m;
    while(sscanf(p," %d-%d%n",&u,&v,&m)==2){ adj[u]|=1u<<v; adj[v]|=1u<<u; p+=m; }
    uint64_t d=count_dom();
    int odd=is_odd();
    int deg[MAXN]; for(int x=0;x<n;x++) deg[x]=__builtin_popcount(adj[x]);
    int U=0; for(int x=0;x<n;x++) if(deg[x]==n-1) U=1; if(n==1) U=1;
    int pend=0,pnu=0,r3=0;
    for(int x=0;x<n;x++){ if(deg[x]==1){ pend++; int y=__builtin_ctz(adj[x]); if(deg[y]<n-1) pnu++; } if(deg[x]==2) r3=1; }
    int fl=1;
    if(!strcmp(set,"W")) fl=!U;
    else if(!strcmp(set,"P")) fl=pend>0;
    else if(!strcmp(set,"S")) fl=pnu>0;
    else if(!strcmp(set,"S2")) fl=pnu>=2;
    else if(!strcmp(set,"S3")) fl=r3;
    else if(!strcmp(set,"G")) fl=(pnu>=2)&&r3;
    if(d==k && !odd && fl) ok++; else { bad++; if(bad<20) printf("BAD %s %s n=%d k=%llu d=%llu odd=%d flags=%d\n",tag,set,n,k,(unsigned long long)d,odd,fl); }
  }
  printf("witnesses ok=%ld bad=%ld\n",ok,bad);
  return bad?1:0;
}
