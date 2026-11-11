/* tower_reach.c - E1 / E2 of the S22 pre-registration: reach of the sign-flip towers.
 * State (r, g, M, E), g = H(B) + r.  Level map (delta in {1,2}, x >= 0, m >= 0), alpha = 2^delta+1, gamma = 3^delta+2^delta:
 *   r+ = -r,  g+ = 2^(x+1+m+delta) M,  M+ = alpha 2^m g + 2 gamma 2^x 3^m M - [x=0] E,  E+ = alpha 2^x M + [m=0](2r - g).
 * usage: tower_reach MODE XMAX MMAXEXP DEPTH CMAX [< states]
 *   MODE 1: pure towers (start states of E1, k <= 10);  MODE 2: start states "r g M E" read from stdin.
 * Output: for each depth d the number of odd c <= CMAX (c = odd part of M) seen at depth <= d with offset +1 and with offset -1;
 * list of (sign, c, least v2(M)) to stdout in the form "S sign c v2 depth". */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef __int128 i128;
static int XMAX, MMAX, DEPTH; static long CMAX;
static unsigned char *seenP, *seenN; static int *v2P, *v2N; static signed char *dP, *dN;
static long pow2M = 0;            /* number of offset -1 states with M a power of two (R3) */
static i128 p3[40];
static void mark(int r, i128 M, int depth){ if(M<=0) return; int v=0; while((M&1)==0){M>>=1;v++;}
  if(r==-1 && M==1) pow2M++;
  if(M>CMAX) return; long c=(long)M; long idx=c>>1;
  if(r==1){ if(!seenP[idx]||v<v2P[idx]){ if(!seenP[idx]) dP[idx]=depth; seenP[idx]=1; v2P[idx]=v; } if(depth<dP[idx]) dP[idx]=depth; }
  else    { if(!seenN[idx]||v<v2N[idx]){ if(!seenN[idx]) dN[idx]=depth; seenN[idx]=1; v2N[idx]=v; } if(depth<dN[idx]) dN[idx]=depth; } }
static void rec(int r, i128 g, i128 M, i128 E, int depth){
  mark(r,M,depth); if(depth==DEPTH) return;
  for(int delta=1;delta<=2;delta++){ int alpha=(1<<delta)+1, gamma=(delta==1)?5:13;
    for(int m=0;m<=MMAX;m++) for(int x=0;x<=XMAX;x++){
      i128 g2=((i128)1<<(x+1+m+delta))*M;
      i128 M2=(i128)alpha*((i128)1<<m)*g + (i128)2*gamma*((i128)1<<x)*p3[m]*M - (x==0?E:0);
      i128 E2=(i128)alpha*((i128)1<<x)*M + (m==0?(2*r-g):0);
      if(depth+1==DEPTH) mark(-r,M2,depth+1); else rec(-r,g2,M2,E2,depth+1); } } }
int main(int argc,char**argv){ int mode=atoi(argv[1]); XMAX=atoi(argv[2]); MMAX=atoi(argv[3]); DEPTH=atoi(argv[4]); CMAX=atol(argv[5]);
  long n=CMAX/2+1; seenP=calloc(n,1); seenN=calloc(n,1); v2P=calloc(n,sizeof(int)); v2N=calloc(n,sizeof(int)); dP=calloc(n,1); dN=calloc(n,1);
  p3[0]=1; for(int i=1;i<40;i++) p3[i]=p3[i-1]*3;
  long starts=0;
  if(mode==1){ for(int k=0;k<=10;k++){ rec(1,2,(i128)1<<k,(k==0),0); starts++; } rec(1,2,0,0,0); starts++; }
  else { long long r,g,M,E; while(scanf("%lld %lld %lld %lld",&r,&g,&M,&E)==4){ rec((int)r,g,M,E,0); starts++; } }
  for(int d=0;d<=DEPTH;d++){ long a=0,b=0; for(long i=0;i<n;i++){ if(seenP[i]&&dP[i]<=d) a++; if(seenN[i]&&dN[i]<=d) b++; }
    fprintf(stderr,"depth <= %d: offset +1: %ld odd c <= %ld ; offset -1: %ld\n",d,a,CMAX,b); }
  fprintf(stderr,"start states %ld ; offset -1 states with M a power of two: %ld\n",starts,pow2M);
  for(long i=0;i<n;i++){ if(seenP[i]) printf("S +1 %ld %d %d\n",2*i+1,v2P[i],dP[i]); if(seenN[i]) printf("S -1 %ld %d %d\n",2*i+1,v2N[i],dN[i]); }
  return 0; }
