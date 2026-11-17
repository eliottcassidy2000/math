/* Wave 59: per class (orders <= 8): n, H, sigma (largest prime-quotient order; 1 = transitive), |Aut|, self-converse */
/* For each class (orders <= 8): H, sigma = largest prime-quotient order in the modular decomposition (1 if transitive),
   |Aut| (brute force over score-respecting permutations), self-converse flag. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static int n; static int out[8];
static u64 Hcount(int m, const int*o){ static u64 dp[256][8]; for(int S=0;S<(1<<m);S++) for(int v=0;v<m;v++) dp[S][v]=0; for(int v=0;v<m;v++) dp[1<<v][v]=1;
  for(int S=1;S<(1<<m);S++) for(int v=0;v<m;v++) if(dp[S][v]) for(int w=0;w<m;w++) if(!(S>>w&1)&&(o[v]>>w&1)) dp[S|1<<w][w]+=dp[S][v];
  u64 h=0; for(int v=0;v<m;v++) h+=dp[(1<<m)-1][v]; return h; }
static int ismodule(int M, int V, const int*o){ /* M subset of V */ for(int x=0;x<8;x++) if((V>>x&1)&&!(M>>x&1)){ int a=-1; for(int y=0;y<8;y++) if(M>>y&1){ int b=o[x]>>y&1; if(a<0) a=b; else if(a!=b) return 0; } } return 1; }
static int strongset(int V, const int*o){ int s=__builtin_ctz(V); int r=1<<s, fr=r; while(fr){ int nf=0; for(int v=0;v<8;v++) if(fr>>v&1) nf|=o[v]&V; nf&=~r; r|=nf; fr=nf; } if(r!=V) return 0;
  r=1<<s; fr=r; while(fr){ int nf=0; for(int v=0;v<8;v++) if(fr>>v&1){ for(int u=0;u<8;u++) if((V>>u&1)&&(o[u]>>v&1)) nf|=1<<u; } nf&=~r; r|=nf; fr=nf; } return r==V; }
static int sigma(int V, const int*o){ int c=__builtin_popcount(V); if(c<=1) return 1;
  if(!strongset(V,o)){ /* split off strong components: find a vertex set closed... use source component: vertices reaching all? simple: recurse on SCCs */
    /* compute SCCs by reachability */ int rem=V, best=1; while(rem){ int s=__builtin_ctz(rem); int f=1<<s, fr=f; while(fr){ int nf=0; for(int v=0;v<8;v++) if(fr>>v&1) nf|=o[v]&V; nf&=~f; f|=nf; fr=nf; }
        int b=1<<s; fr=b; while(fr){ int nf=0; for(int v=0;v<8;v++) if(fr>>v&1){ for(int u=0;u<8;u++) if((V>>u&1)&&(o[u]>>v&1)) nf|=1<<u; } nf&=~b; b|=nf; fr=nf; }
        int scc=f&b; int sg=sigma(scc,o); if(sg>best) best=sg; rem&=~scc; } return best; }
  /* strong: maximal proper modules partition V */
  int mods[256], nm=0; for(int M=V; M; M=(M-1)&V){ if(M==V||__builtin_popcount(M)<1) continue; if(ismodule(M,V,o)) mods[nm++]=M; }
  int parts[8], np=0, cover=0;
  for(int i=0;i<nm;i++){ int M=mods[i]; int maximal=1; for(int j=0;j<nm;j++) if(mods[j]!=M && (mods[j]&M)==M){ maximal=0; break; } if(maximal && !(cover&M)){ parts[np++]=M; cover|=M; } }
  int best=np; for(int i=0;i<np;i++){ int sg=sigma(parts[i],o); if(sg>best) best=sg; } return best; }
static int perm[8], used[8], cnt_aut, conv_found; static int score[8];
static void autrec(int pos, int mode){ /* mode 0: automorphisms; 1: anti-automorphisms (converse isomorphisms) */
  if(pos==n){ if(mode==0) cnt_aut++; else conv_found=1; return; }
  if(mode==1 && conv_found) return;
  for(int v=0;v<n;v++){ if(used[v]) continue; if(mode==0 && score[v]!=score[pos]) continue; if(mode==1 && score[v]!=n-1-score[pos]) continue;
    int ok=1; for(int q=0;q<pos&&ok;q++){ int a=out[pos]>>q&1, b=out[v]>>perm[q]&1; if(mode==0){ if(a!=b) ok=0; } else { if(a==b) ok=0; } }
    if(!ok) continue; used[v]=1; perm[pos]=v; autrec(pos+1,mode); used[v]=0; } }
int main(void){ FILE*fi=fopen("w59_reps8.txt","r"); char line[256];
  while(fgets(line,sizeof line,fi)){ char*p=line; n=strtol(p,&p,10); for(int i=0;i<n;i++) out[i]=strtol(p,&p,10);
    for(int i=0;i<n;i++) score[i]=__builtin_popcount(out[i]);
    u64 h=Hcount(n,out); int sg=sigma((1<<n)-1,out);
    memset(used,0,sizeof used); cnt_aut=0; autrec(0,0); memset(used,0,sizeof used); conv_found=0; autrec(0,1);
    printf("%d %llu %d %d %d\n",n,h,sg,cnt_aut,conv_found); }
  return 0; }
