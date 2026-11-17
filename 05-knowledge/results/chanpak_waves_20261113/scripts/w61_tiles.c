/* Wave 61. */
/* Signed tiles (signed_ordinal_unit_centers_20261117 §1): tile D with port masks Ix, Ia (vertices beating x, a).
   Endpoint match: the start and end vectors of Hamiltonian paths and the (end of first, start of second) vectors of
   ordered spanning two-path covers have equal sums on Ix and Ia.  Charge 2s = Ix^T (M - M^T) Ia.
   Census over all classes of orders 7 and 8: matched charged pairs (s != 0), split by 3 | H(D); unit charges s = +-1. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
static ll endE[256][8], startS[256][8], Hs[256];
int main(int argc,char**argv){ FILE*fi=fopen("w59_reps8.txt","r"); char line[256];
  long cnt[9][2]={{0}}, unit[9][2]={{0}}, classes_with[9][2]={{0}}; ll maxs[9]={0};
  long unitcls7=0; int printed=0; long bez[9]={0}, unitC3[9]={0}, unitnoC3[9]={0};
  while(fgets(line,sizeof line,fi)){ char*p=line; int n=strtol(p,&p,10); int o[8]; for(int i=0;i<n;i++) o[i]=strtol(p,&p,10); if(n<7) continue;
    int N=1<<n; for(int S=0;S<N;S++) for(int v=0;v<n;v++){ endE[S][v]=0; startS[S][v]=0; }
    for(int v=0;v<n;v++){ endE[1<<v][v]=1; startS[1<<v][v]=1; }
    for(int S=1;S<N;S++) for(int v=0;v<n;v++){ if(!(S>>v&1)) continue;
      if(endE[S][v]) for(int w=0;w<n;w++) if(!(S>>w&1)&&(o[v]>>w&1)) endE[S|1<<w][w]+=endE[S][v];
      if(startS[S][v]) for(int w=0;w<n;w++) if(!(S>>w&1)&&(o[w]>>v&1)) startS[S|1<<w][w]+=startS[S][v]; }
    ll H=0; for(int v=0;v<n;v++) H+=endE[N-1][v];
    ll st[8],en[8],e1[8]={0},s2[8]={0},M[8][8]; for(int v=0;v<n;v++){ st[v]=startS[N-1][v]; en[v]=endE[N-1][v]; }
    for(int u=0;u<n;u++) for(int v=0;v<n;v++){ M[u][v]=0; if(u==v) continue; for(int S=1;S<N;S++) if((S>>u&1)&&!(S>>v&1)) M[u][v]+=endE[S][u]*startS[(N-1)&~S][v]; }
    for(int u=0;u<n;u++) for(int v=0;v<n;v++){ e1[u]+=M[u][v]; s2[v]+=M[u][v]; }
    static ll sst[256],sen[256],se1[256],ss2[256]; static ll A[8][256];
    for(int m=0;m<N;m++){ sst[m]=sen[m]=se1[m]=ss2[m]=0; for(int v=0;v<n;v++) if(m>>v&1){ sst[m]+=st[v]; sen[m]+=en[v]; se1[m]+=e1[v]; ss2[m]+=s2[v]; } }
    for(int u=0;u<n;u++) for(int m=0;m<N;m++){ ll t=0; for(int v=0;v<n;v++) if(m>>v&1) t+=M[u][v]-M[v][u]; A[u][m]=t; }
    int d3=(H%3==0); long here=0, hereunit=0;
    for(int ix=0;ix<N;ix++) for(int ia=0;ia<N;ia++){
      if(sst[ix]!=sst[ia]||sen[ix]!=sen[ia]||se1[ix]!=se1[ia]||ss2[ix]!=ss2[ia]) continue;
      ll two_s=0; for(int u=0;u<n;u++) if(ix>>u&1) two_s+=A[u][ia];
      if(two_s==0) continue; if(two_s%2){ printf("odd 2s?!\n"); return 1; }
      ll s=two_s/2; cnt[n][d3]++; here++; if(s==1||s==-1){ unit[n][d3]++; hereunit++; }
      if((s>0?s:-s)>maxs[n]) maxs[n]=(s>0?s:-s);
      if(!d3 && (3*s-H==1 || 3*s-H==-1 || 3*s+H==1 || 3*s+H==-1)){ bez[n]++; printf("  Bezout-3 companion n=%d H=%lld s=%lld (3s-H=%lld) rows",n,H,s,3*s-H); for(int i=0;i<n;i++) printf(" %d",o[i]); printf(" Ix=%d Ia=%d\n",ix,ia); }
      if(!d3 && (s==1||s==-1)){ int hasC3=0; for(int a=0;a<n&&!hasC3;a++) for(int b=0;b<n&&!hasC3;b++) for(int c=0;c<n&&!hasC3;c++){ if(a==b||b==c||a==c) continue; if(!((o[a]>>b&1)&&(o[b]>>c&1)&&(o[c]>>a&1))) continue; int m=(1<<a)|(1<<b)|(1<<c); int ok=1; for(int x=0;x<n&&ok;x++) if(!(m>>x&1)){ int ab=o[x]>>a&1, bb=o[x]>>b&1, cb=o[x]>>c&1; if(ab!=bb||bb!=cb) ok=0; } if(ok) hasC3=1; } if(hasC3) unitC3[n]++; else unitnoC3[n]++; }
      if(n==7 && !d3 && (s==1||s==-1) && printed<6){ printed++; printf("  unit tile n=7 H=%lld s=%lld rows",H,s); for(int i=0;i<n;i++) printf(" %d",o[i]); printf(" Ix=%d Ia=%d\n",ix,ia); } }
    if(here) classes_with[n][d3]++; }
  for(int n=7;n<=8;n++) printf("order %d: matched charged mask pairs with 3 | H: %ld (classes %ld), with 3 not dividing H: %ld (classes %ld); unit charges: %ld / %ld; max |s| %lld\n",
     n,cnt[n][1],classes_with[n][1],cnt[n][0],classes_with[n][0],unit[n][1],unit[n][0],maxs[n]);
  for(int n=7;n<=8;n++) printf("order %d, 3 not dividing H: Bezout-3 pairs (3s - h = +-1 up to sign) %ld; unit pairs with an internal C3 module %ld, without %ld\n",n,bez[n],unitC3[n],unitnoC3[n]);
  return 0; }
