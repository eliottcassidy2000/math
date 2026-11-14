// Oracle study: which statistic of a prefix subtree predicts that it contains a word with H = N?
// Enumerate all two-hub words of length L (both e); group leaves by their first K letters.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define NS 36
typedef unsigned long long u64;
static int idx(int o,int a,int b,int d){return ((o*3+a)*3+b)*2+d;}
static int T[4][NS][3], Tn[4][NS], fin[2][NS];
static void build(void){
  for(int let=0;let<4;let++){ int xi=let>>1, yi=let&1;
    for(int o=0;o<2;o++)for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int d=0;d<2;d++){ int s=idx(o,a,b,d), k=0;
      int h1=o?yi:xi, h2=o?xi:yi;
      T[let][s][k++]=idx(o,h1?2:1,b,d);
      if(b||h1) T[let][s][k++]=idx(o,a,h2?2:1,d);
      if(d||h2) T[let][s][k++]=idx(o,a,b,1);
      Tn[let][s]=k; } }
  for(int e=0;e<2;e++) for(int o=0;o<2;o++)for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int d=0;d<2;d++){
      int h12=o?!e:e; fin[e][idx(o,a,b,d)] = !(a==2||b==2||(b==0&&!h12)); } }
static int L,K; static u64 *Hs; static int *Pf; static long nleaf=0;
static void dfs(int k,u64*v,int pref){
  if(k==L){ for(int e=0;e<2;e++){ u64 H=0; for(int s=0;s<NS;s++) if(fin[e][s]) H+=v[s]; Hs[nleaf]=H; Pf[nleaf]=pref; nleaf++; } return; }
  for(int let=0;let<4;let++){ u64 w[NS]; memset(w,0,sizeof w);
    for(int s=0;s<NS;s++) if(v[s]) for(int t=0;t<Tn[let][s];t++) w[T[let][s][t]]+=v[s];
    dfs(k+1,w,k<K?pref*4+let:pref); }
}
static int cmpu(const void*a,const void*b){ u64 x=*(u64*)a,y=*(u64*)b; return x<y?-1:x>y; }
int main(int argc,char**argv){
  L=atoi(argv[1]); K=atoi(argv[2]); build();
  long n=2L<<(2*L); Hs=malloc(sizeof(u64)*n); Pf=malloc(sizeof(int)*n);
  u64 v[NS]; memset(v,0,sizeof v); v[idx(0,0,0,0)]=1; v[idx(1,0,0,0)]=1; dfs(0,v,0);
  int P=1<<(2*K);
  // per-prefix sorted H lists
  long *cnt=calloc(P,sizeof(long)); for(long i=0;i<n;i++) cnt[Pf[i]]++;
  u64 **lst=malloc(sizeof(u64*)*P); long *fill=calloc(P,sizeof(long));
  for(int p=0;p<P;p++) lst[p]=malloc(sizeof(u64)*(cnt[p]+1));
  for(long i=0;i<n;i++) lst[Pf[i]][fill[Pf[i]]++]=Hs[i];
  double *mu=calloc(P,8),*sd=calloc(P,8),*lmu=calloc(P,8),*lsd=calloc(P,8); u64 *mn=calloc(P,8),*mx=calloc(P,8);
  int MODS[]={3,4,5,8,9,16,27}; int NM=7; double *res[7];
  for(int j=0;j<NM;j++) res[j]=calloc((size_t)P*MODS[j],8);
  for(int p=0;p<P;p++){ qsort(lst[p],cnt[p],8,cmpu); double s=0,s2=0,l=0,l2=0; mn[p]=lst[p][0]; mx[p]=lst[p][cnt[p]-1];
    for(long i=0;i<cnt[p];i++){ double h=(double)lst[p][i]; s+=h; s2+=h*h; double lh=log(h>0?h:0.5); l+=lh; l2+=lh*lh;
      for(int j=0;j<NM;j++) res[j][(size_t)p*MODS[j]+lst[p][i]%MODS[j]]+=1; }
    mu[p]=s/cnt[p]; sd[p]=sqrt(fmax(s2/cnt[p]-mu[p]*mu[p],1e-9)); lmu[p]=l/cnt[p]; lsd[p]=sqrt(fmax(l2/cnt[p]-lmu[p]*lmu[p],1e-9));
    for(int j=0;j<NM;j++) for(int c=0;c<MODS[j];c++) res[j][(size_t)p*MODS[j]+c]/=cnt[p]; }
  // horizon: first odd miss above 21
  u64 *all=malloc(sizeof(u64)*n); memcpy(all,Hs,sizeof(u64)*n); qsort(all,n,8,cmpu);
  u64 X=23; { long i=0; while(1){ while(i<n && all[i]<X) i++; if(i<n && all[i]==X) X+=2; else break; } }
  // evaluate heuristics over odd N in [lo, X)
  u64 lo=X/9; if(lo<23) lo=23; lo|=1;
  const char *names[]={"random","logcenter(min,max)","lognormal density","normal density","residue mod 9 (oracle)","residue mod 27 (oracle)","lognormal x res9","lognormal x res27","oracle best"};
  int NH=9; long succ[9]={0}; long tot=0;
  for(u64 N=lo; N<X; N+=2){ if(N==21) continue; tot++;
    double sc[9][1024]; long sol[1024];
    for(int p=0;p<P;p++){
      // solutions in prefix p
      u64 *a=lst[p]; long lo2=0,hi2=cnt[p]; while(lo2<hi2){ long md=(lo2+hi2)/2; if(a[md]<N) lo2=md+1; else hi2=md; }
      long c=0; while(lo2+c<cnt[p] && a[lo2+c]==N) c++; sol[p]=c;
      int feas=(N>=mn[p] && N<=mx[p]);
      double z=(log((double)N)-lmu[p])/lsd[p]; double ln=feas? cnt[p]*exp(-0.5*z*z)/(lsd[p]*N):0;
      double zz=((double)N-mu[p])/sd[p]; double nd=feas? cnt[p]*exp(-0.5*zz*zz)/sd[p]:0;
      double r9=feas? res[4][(size_t)p*9+N%9]:0, r27=feas? res[6][(size_t)p*27+N%27]:0;
      double lc=feas? -fabs(log((double)N)-0.5*(log((double)mn[p])+log((double)mx[p]))):-1e300;
      sc[0][p]=feas? (double)((p*2654435761u+N*40503u)%1000003):-1; sc[1][p]=lc; sc[2][p]=ln; sc[3][p]=nd; sc[4][p]=r9*cnt[p]; sc[5][p]=r27*cnt[p];
      sc[6][p]=ln*r9; sc[7][p]=ln*r27; sc[8][p]=(double)sol[p];
    }
    for(int h=0;h<NH;h++){ int best=-1; double bv=-1e301; for(int p=0;p<P;p++) if(sc[h][p]>bv){bv=sc[h][p];best=p;} if(best>=0 && sol[best]>0) succ[h]++; }
  }
  printf("L=%d, K=%d prefix letters (%d prefixes), horizon X=%llu, targets odd N in [%llu,%llu): %ld\n",L,K,P,X,lo,X,tot);
  for(int h=0;h<NH;h++) printf("  %-26s top choice has a solution: %.1f%%\n",names[h],100.0*succ[h]/tot);
  return 0;
}
