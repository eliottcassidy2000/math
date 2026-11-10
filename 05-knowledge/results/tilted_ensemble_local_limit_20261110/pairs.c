// pairs.c - two-dimensional census of deletion pairs (H(T), H(T+z)) over all tournaments T of order n
// and all 2^n out-sets of the new vertex z.  Records, for each h < HM, the bitset of k = (H(T+z)-h)/2 < KM.
// Also the Z-subset (out-degree one extensions) separately.
// Output: binary dump of two bitset tables to files given on the command line (OR-merged if they exist).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64; typedef long long i64;
#define MAXN 11
#define HM 1024      /* h < HM (odd) */
#define KM 4096      /* k = (n-h)/2 < KM */
static int n;
static u64 f[1<<MAXN][MAXN], g[1<<MAXN][MAXN], hS[1<<MAXN];
static int outm[MAXN], inm[MAXN];
static unsigned char *tabA, *tabZ;   /* [HM/2][KM/8] */
#define SETB(t,h,k) (t[((size_t)((h)>>1))*(KM/8) + ((k)>>3)] |= (unsigned char)(1u<<((k)&7)))
int main(int argc, char **argv){
  n = atoi(argv[1]); const char *fa = argv[2], *fz = argv[3];
  size_t sz = (size_t)(HM/2)*(KM/8);
  tabA = calloc(sz,1); tabZ = calloc(sz,1);
  FILE *fp;
  if ((fp=fopen(fa,"rb"))){ if(fread(tabA,1,sz,fp)!=sz){} fclose(fp);} 
  if ((fp=fopen(fz,"rb"))){ if(fread(tabZ,1,sz,fp)!=sz){} fclose(fp);} 
  char line[256]; int full=(1<<n)-1; u64 nt=0;
  static u64 M[MAXN][MAXN];
  while (fgets(line,sizeof line,stdin)){
    if (line[0]!='0'&&line[0]!='1') continue;
    int k=0; for(int i=0;i<n;i++){outm[i]=0;inm[i]=0;}
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){ if(line[k++]=='1'){outm[i]|=1<<j;inm[j]|=1<<i;} else {outm[j]|=1<<i;inm[i]|=1<<j;} }
    nt++; hS[0]=1;
    for(int S=1;S<=full;S++){ u64 tot=0;
      for(int v=0;v<n;v++){ if(!(S>>v&1)){f[S][v]=g[S][v]=0;continue;}
        int R=S&~(1<<v); if(R==0){f[S][v]=g[S][v]=1;tot+=1;continue;}
        u64 a=0,b=0; int m=R&inm[v]; while(m){int w=__builtin_ctz(m);m&=m-1;a+=f[R][w];}
        m=R&outm[v]; while(m){int w=__builtin_ctz(m);m&=m-1;b+=g[R][w];}
        f[S][v]=a;g[S][v]=b;tot+=a; }
      hS[S]=tot; }
    u64 H=hS[full]; if (H>=HM) continue;
    for(int j=0;j<n;j++)for(int i=0;i<n;i++){ if(i==j){M[j][i]=0;continue;}
      u64 a=0; int rest=full&~(1<<i)&~(1<<j);
      for(int S0=rest;;S0=(S0-1)&rest){int S=S0|(1<<j); a+=f[S][j]*g[full^S][i]; if(S0==0)break;}
      M[j][i]=a; }
    for(int X=0;X<=full;X++){ u64 val=0;
      for(int i=0;i<n;i++){ if(X>>i&1) val+=g[full][i]; else { val+=f[full][i]; int m=X; while(m){int t=__builtin_ctz(m);m&=m-1;val+=M[i][t];} } }
      u64 kk=(val-H)/2; if(kk<KM){ SETB(tabA,H,kk); if(__builtin_popcount(X)==1) SETB(tabZ,H,kk); }
    }
  }
  fp=fopen(fa,"wb"); fwrite(tabA,1,sz,fp); fclose(fp);
  fp=fopen(fz,"wb"); fwrite(tabZ,1,sz,fp); fclose(fp);
  fprintf(stderr,"n=%d tournaments=%llu\n",n,nt); return 0;
}
