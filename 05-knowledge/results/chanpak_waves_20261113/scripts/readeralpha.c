// One-hub source U followed by a chain of transitive readers (c a^i b^j) drawn from a given alphabet.
// Notation: these are the repository's O/Z/J root-change words (THM-4904, THM-5148, HYP-9494): a = O, b = Z,
// c = J.  The reader blocks c a^i b^j = J O^i Z^j are not transitive; their counts are (2^i-1)(2^j-1)+2.
// usage: readeralpha X i1 j1 i2 j2 ...
#include <stdio.h>
#include <stdlib.h>
typedef long long i64;
static i64 X; static unsigned char *cov; static int NR; static int RI[256],RJ[256]; static i64 M[256][3][3];
static i64 Hc(i64 U){ if(U==0) return 1; int k=64-__builtin_clzll((unsigned long long)U); i64 s=1;
  for(int t=0;t<k;t++){ int bt=(U>>t)&1; if(!bt) s+=U>>(t+1); if(t+1<k && ((U>>(t+1))&1) && !bt) s+=1; } return s; }
static void st_a(i64*s){ s[1]=0; s[2]=2*s[2]; }
static void st_b(i64*s){ i64 H=s[0],E=s[1],L=s[2],h=s[3]; s[0]=H+L-E; s[1]=h; s[2]=2*L-h; }
static void st_c(i64*s){ i64 H=s[0],E=s[1],L=s[2]; s[0]=H+L-E; s[1]=H-E; s[2]=H+L; s[3]=H; }
static void dfs(i64 H,i64 E,i64 L,int depth){
  if(H<=X) cov[H>>1]=1;
  if(depth>60) return;
  for(int r=0;r<NR;r++){
    i64 nH=M[r][0][0]*H+M[r][0][1]*E+M[r][0][2]*L; if(nH>X || nH<=0) continue;
    i64 nE=M[r][1][0]*H+M[r][1][1]*E+M[r][1][2]*L, nL=M[r][2][0]*H+M[r][2][1]*E+M[r][2][2]*L;
    if(nH==H && nL<=L) continue;
    dfs(nH,nE,nL,depth+1);
  }
}
int main(int argc,char**argv){
  X=atoll(argv[1]); NR=(argc-2)/2; for(int r=0;r<NR;r++){ RI[r]=atoi(argv[2+2*r]); RJ[r]=atoi(argv[3+2*r]); }
  for(int r=0;r<NR;r++) for(int col=0;col<3;col++){ i64 s[4]={0,0,0,0}; s[col]=1; st_c(s); for(int t=0;t<RI[r];t++) st_a(s); for(int t=0;t<RJ[r];t++) st_b(s); for(int row=0;row<3;row++) M[r][row][col]=s[row]; }
  cov=calloc(X/2+2,1);
  for(i64 U=0;U<=2*X;U++){ i64 H=Hc(U); if(H<=X) cov[H>>1]=1; if(U>X) continue; i64 E=(U&1)?0:1; dfs(H,E,1+U,0); }
  long miss=0,last=0; int shown=0; printf("X=%lld alphabet:",X); for(int r=0;r<NR;r++) printf(" (%d,%d)",RI[r],RJ[r]); printf("\n misses:");
  for(i64 n=1;n<=X;n+=2) if(!cov[n>>1]){ miss++; last=n; if(shown<30){printf(" %lld",n);shown++;} }
  printf("\n total %ld, last %ld\n",miss,last); return 0; }
