// allpairs.c: for every tournament of order n on stdin (gentourng upper-triangle format) and EVERY vertex a,
// compute h = h(T), e = e_a(T), B = B_a(T) = sum_{A subset V\a} h(T[A]) s_a(T[V\A]); print "h e B strong" per (T,a).
// (The one-exit extension T + x (x -> a) has h - e + B Hamiltonian paths; the C2 child has h + 2B.)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static int n; static u64 adj[16], radj[16];
static u64 *f, *g; // f[S*n+v]: # ham paths of T[S] ending at v; g[S*n+v]: # ham paths of T[S] starting at v
static u64 hS[1<<14];
static void dp(u64 *tab, u64 *A){ u64 N=1ULL<<n; memset(tab,0,sizeof(u64)*N*n);
  for(int v=0;v<n;v++) tab[(1ULL<<v)*n+v]=1;
  for(u64 S=1;S<N;S++) for(int v=0;v<n;v++){ u64 c=tab[S*n+v]; if(!c) continue; u64 nx=A[v]&~S; while(nx){int u=__builtin_ctzll(nx); nx&=nx-1; tab[(S|(1ULL<<u))*n+u]+=c;} } }
static int is_strong(void){ u64 all=(1ULL<<n)-1, seen=1, fr=1; while(fr){u64 nf=0; for(int v=0;v<n;v++) if(fr>>v&1) nf|=adj[v]; nf&=~seen; seen|=nf; fr=nf;} if(seen!=all) return 0; seen=1; fr=1; while(fr){u64 nf=0; for(int v=0;v<n;v++) if(fr>>v&1) nf|=radj[v]; nf&=~seen; seen|=nf; fr=nf;} return seen==all; }
int main(int argc,char**argv){ n=atoi(argv[1]); u64 N=1ULL<<n; f=malloc(sizeof(u64)*N*n); g=malloc(sizeof(u64)*N*n);
  char line[4096]; u64 total=0;
  while(fgets(line,sizeof line,stdin)){ int L=strlen(line); while(L&&(line[L-1]=='\n'||line[L-1]=='\r')) line[--L]=0; if(L!=n*(n-1)/2) continue;
    memset(adj,0,sizeof adj); memset(radj,0,sizeof radj); int idx=0;
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if(line[idx++]=='1'){adj[i]|=1ULL<<j; radj[j]|=1ULL<<i;} else {adj[j]|=1ULL<<i; radj[i]|=1ULL<<j;} }
    total++; int st=is_strong();
    dp(f,adj); dp(g,radj); // g computed on reversed graph: g[S][v] = # ham paths of T[S] ending at v in reversed = starting at v in T
    for(u64 S=1;S<N;S++){ u64 s=0; for(int v=0;v<n;v++) s+=f[S*n+v]; hS[S]=s; } hS[0]=1;
    u64 h=hS[N-1];
    for(int a=0;a<n;a++){ u64 e=f[(N-1)*n+a]; u64 B=0; u64 rest=(N-1)&~(1ULL<<a);
      // enumerate A subset of rest
      for(u64 A=rest;;A=(A-1)&rest){ u64 W=(N-1)&~A; B+=hS[A]*g[W*n+a]; if(A==0) break; }
      printf("%llu %llu %llu %d\n",h,e,B,st); }
  }
  fprintf(stderr,"order %d: %llu tournaments\n",n,total); return 0; }
