// Count representations of N in a root-change (O/Z/J) sub-language (letters a=O, b=Z, c=J):
// one-apex source U (any length) + chain of readers (J O^i Z^j), i,j <= MAXIJ, at most MAXR readers.
// Representation = (chain, U) with U canonical (U even or U==0: trailing O's are free), counted exactly.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef __int128 i128; typedef unsigned long long u64;
static i128 Hc(i128 U){ if(U==0) return 1; int k=0; for(i128 t=U;t;t>>=1) k++; i128 s=1;
  for(int t=0;t<k;t++){ int bt=(int)((U>>t)&1); if(!bt) s+=U>>(t+1); if(t+1<k && ((U>>(t+1))&1) && !bt) s+=1; } return s; }
static void st_a(i128*s){ s[1]=0; s[2]=2*s[2]; }
static void st_b(i128*s){ i128 H=s[0],E=s[1],L=s[2],h=s[3]; s[0]=H+L-E; s[1]=h; s[2]=2*L-h; }
static void st_c(i128*s){ i128 H=s[0],E=s[1],L=s[2]; s[0]=H+L-E; s[1]=H-E; s[2]=H+L; s[3]=H; }
static i128 M[16][16][3][3]; static int MAXR,MAXIJ; static i128 N0; static u64 nchains_ok;
// brute count over U directly for small U ranges (exact, no recursion subtleties): enumerate U in [lo,hi]
static u64 count_direct(i128 A,i128 C,i128 B){
  // N = A H(U) + C E(U) + B (1+U); E(U)=1-(U&1). Canonical: count all U (each U is a distinct tournament word).
  u64 c=0; if(A+B+C>N0) return 0;
  i128 lo,hi;
  if(B>0){ lo=(N0-A-B-(C>0?C:0))/(A+B)-2; if(lo<0) lo=0; hi=(N0-A-B-(C<0?C:0))/B+2; }
  else { lo=0; hi=2*N0+4; }
  if(hi-lo>4000000) return (u64)-1;
  for(i128 U=lo;U<=hi;U++){ i128 E=(U&1)?0:1; if(A*Hc(U)+C*E+B*(1+U)==N0) c++; }
  return c;
}
static u64 total; static int overflow;
static void chains(i128 A,i128 C,i128 B,int r){
  if(A+B+C>N0) return;
  if(r>0 || 1){ u64 c=count_direct(A,C,B); if(c==(u64)-1) overflow=1; else { total+=c; if(c) nchains_ok++; } }
  if(r>=MAXR) return;
  for(int i=0;i<=MAXIJ;i++) for(int j=0;j<=MAXIJ;j++){ i128 (*R)[3]=M[i][j];
    i128 nA=A*R[0][0]+C*R[1][0]+B*R[2][0], nC=A*R[0][1]+C*R[1][1]+B*R[2][1], nB=A*R[0][2]+C*R[1][2]+B*R[2][2];
    chains(nA,nC,nB,r+1); }
}
int main(int argc,char**argv){
  i128 LO=atoll(argv[1]); int cnt=atoi(argv[2]); MAXR=atoi(argv[3]); MAXIJ=atoi(argv[4]);
  for(int i=0;i<=MAXIJ;i++) for(int j=0;j<=MAXIJ;j++) for(int col=0;col<3;col++){ i128 s[4]={0,0,0,0}; s[col]=1; st_c(s); for(int t=0;t<i;t++) st_a(s); for(int t=0;t<j;t++) st_b(s); for(int row=0;row<3;row++) M[i][j][row][col]=s[row]; }
  u64 mn=(u64)-1, mx=0; double sum=0; int n=0;
  for(int t=0;t<cnt;t++){ N0=LO+2*t+1-(LO&1?1:0); if(!(N0&1)) N0++; total=0; nchains_ok=0; overflow=0;
    chains(1,0,0,0);
    // the no-reader case (pure one-apex) has B=0: count_direct scans U up to 2N -- skip if too big
    if(overflow){ printf("N=%lld overflow (no-reader scan too big); counting readers only\n",(long long)N0); }
    if(total<mn) mn=total; if(total>mx) mx=total; sum+=total; n++; }
  printf("N in [%lld..] %d targets, MAXR=%d MAXIJ=%d: representations min %llu mean %.1f max %llu\n",(long long)LO,n,MAXR,MAXIJ,mn,sum/n,mx);
  return 0; }
