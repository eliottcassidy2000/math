// catdec: residue-steered, target-first decoder for caterpillar tournaments (shape 3).
// Notation: these are the repository's O/Z/J root-change words (THM-4904, THM-5148, HYP-9494): a = O, b = Z,
// c = J.  The reader blocks c a^i b^j = J O^i Z^j are not transitive; their counts are (2^i-1)(2^j-1)+2.
//
// A caterpillar word is  x  (c a^i1 b^j1) (c a^i2 b^j2) ... (c a^ir b^jr)
// x = one-hub source word (a = no back arc, b = leaf), each (c a^i b^j) a transitive reader.
// Root-state recurrence (H,E,L,h): a:(H,0,2L,h)  b:(H+L-E,h,2L-h,h)  c:(H+L-E,H-E,H+L,H), start (1,1,1,1).
//
// Decoding N:
//  1. Peel readers from the END: the count is a linear functional phi=(A,C,B) of the source triple
//     (H(U), E(U), 1+U), U = binary value of x (a=1,b=0, MSB first).  N = A H(U) + C E(U) + B (1+U).
//  2. Put the factors into the source: solve A*H(U) + B*U = T by peeling the bits of U from the bottom,
//     using the split law H(2U'+w) = H(U') + [w=0](U' + [U' odd]).  Each peeled bit updates
//     (B,T) -> (2B, T-B) for w=1 or (2B+A, T - A[U' odd]) for w=0: the reader's factor B is carried
//     inside the source problem, and the next bit's parity is a residue requirement.
//  3. Final absolute equality: the recursion closes only at U=0 with T == A exactly (or by direct check
//     once the interval [(T-A)/(A+B), (T-A)/B] holds few integers).  Every success is re-verified by
//     running the 4D recurrence on the produced word.
//
// usage: catdec N [maxreaders] [maxij]     (N up to ~1e30)   or   catdec -r LO HI count [maxreaders] [maxij]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
typedef __int128 i128; typedef unsigned long long u64;

static i128 Hc(i128 U){ // closed form one-hub count
  if(U==0) return 1;
  int k=0; for(i128 t=U;t;t>>=1) k++;
  i128 s=1;
  for(int t=0;t<k;t++){
    int bt=(int)((U>>t)&1);
    if(!bt) s+=U>>(t+1);
    if(t+1<k && ((U>>(t+1))&1) && !bt) s+=1;
  }
  return s;
}
static void p128(i128 x,char*buf){ char t[64]; int n=0; int neg=x<0; if(neg) x=-x; if(x==0) t[n++]='0'; while(x){t[n++]='0'+(int)(x%10); x/=10;} int m=0; if(neg) buf[m++]='-'; while(n) buf[m++]=t[--n]; buf[m]=0; }
static i128 parse128(const char*s){ i128 x=0; while(*s){ x=x*10+(*s-'0'); s++;} return x; }

static u64 nodes; static long long NODECAP;
// solve A*H(U) + B*U == T with parity requirement req (0 any, 1 even, 2 odd). Returns U or -1.
static i128 solveU(i128 A,i128 B,i128 T,int req,int depth){
  nodes++; if(nodes>(u64)NODECAP) return -1;
  if(T<A) return -1;
  if(T==A && req!=2) return 0;                       // U = 0
  if(B>0){
    // A + B*U <= T <= A*(U+1) + B*U  for U >= 1
    i128 lo=(T-A)/(A+B); if(lo<1) lo=1; i128 hi=(T-A)/B;
    if(lo>hi) return -1;
    if(hi-lo<=8){
      for(i128 U=lo;U<=hi;U++){
        if(req==1 && (U&1)) continue; if(req==2 && !(U&1)) continue;
        if(A*Hc(U)+B*U==T) return U;
      }
      return -1;
    }
  }
  if(depth>125) return -1;
  // w = 1 (U odd): U = 2U'+1, H(U)=H(U')  -> A H(U') + 2B U' = T - B
  if(req!=1 && B>0){
    i128 r=solveU(A,2*B,T-B,0,depth+1);
    if(r>=0) return 2*r+1;
  }
  // w = 0 (U even, U>=2): U = 2U', H(U) = H(U') + U' + [U' odd]
  if(req!=2){
    i128 r=solveU(A,2*B+A,T-A,2,depth+1); if(r>=0) return 2*r;      // U' odd
    r=solveU(A,2*B+A,T,1,depth+1);   if(r>=1) return 2*r;             // U' even, U' >= 1 (U'=0 gives U=0)
  }
  return -1;
}

// transitive reader (i,j) acting on the triple (H,E,L)
typedef struct { i128 m[3][3]; } M3;
static void st_a(i128*s){ s[1]=0; s[2]=2*s[2]; }
static void st_b(i128*s){ i128 H=s[0],E=s[1],L=s[2],h=s[3]; s[0]=H+L-E; s[1]=h; s[2]=2*L-h; }
static void st_c(i128*s){ i128 H=s[0],E=s[1],L=s[2]; s[0]=H+L-E; s[1]=H-E; s[2]=H+L; s[3]=H; }
static M3 reader(int i,int j){
  M3 M; for(int col=0;col<3;col++){ i128 s[4]={0,0,0,0}; s[col]=1; st_c(s); for(int t=0;t<i;t++) st_a(s); for(int t=0;t<j;t++) st_b(s);
    for(int row=0;row<3;row++) M.m[row][col]=s[row]; }
  return M;
}
static int MAXR, MAXIJ; static M3 RD[64][64];
static int seqi[32], seqj[32]; static int found; static i128 foundU; static int foundr;
static i128 N0;
// try the source with functional phi=(A,C,B):  N = A H(U) + C E(U) + B (1+U), E(U) = 1 - (U&1)
static long long CAP=40000000LL; static u64 totnodes;
static int try_source(i128 A,i128 C,i128 B){
  if(A<=0||B<0) return 0;
  double est=__builtin_sqrt((double)N0*(double)A)/(double)(B>0?B:1);
  nodes=0; NODECAP=(long long)(6*est)+512; if(NODECAP>CAP) NODECAP=CAP;
  i128 U; int ok=0;
  U=solveU(A,B,N0-B,2,0); if(U>=0){ foundU=U; ok=1; }               // U odd: E=0
  if(!ok){ U=solveU(A,B,N0-B-C,1,0); if(U>=0){ foundU=U; ok=1; } }  // U even: E=1
  totnodes+=nodes; return ok;
}
// enumerate reader chains of total size exactly S (size of (c a^i b^j) = 1+i+j), peeled from the end
static void chains(i128 A,i128 C,i128 B,int r,int S){
  if(found) return;
  if(A+B+C> N0) return;
  if(S==0){ if(try_source(A,C,B)){ found=1; foundr=r; } return; }
  if(r>=MAXR) return;
  for(int sz=1; sz<=S && !found; sz++) for(int i=0;i<sz && !found;i++){ int j=sz-1-i; if(i>MAXIJ||j>MAXIJ) continue;
      M3*M=&RD[i][j];
      i128 nA=A*M->m[0][0]+C*M->m[1][0]+B*M->m[2][0];
      i128 nC=A*M->m[0][1]+C*M->m[1][1]+B*M->m[2][1];
      i128 nB=A*M->m[0][2]+C*M->m[1][2]+B*M->m[2][2];
      seqi[r]=i; seqj[r]=j;
      chains(nA,nC,nB,r+1,S-sz);
  }
}
// verification by the 4D recurrence on the produced word
static i128 verify(i128 U,int r){
  i128 s[4]={1,1,1,1}; int k=0; for(i128 t=U;t;t>>=1) k++;
  for(int t=k-1;t>=0;t--){ if((U>>t)&1) st_a(s); else st_b(s); }
  for(int q=r-1;q>=0;q--){ st_c(s); for(int t=0;t<seqi[q];t++) st_a(s); for(int t=0;t<seqj[q];t++) st_b(s); }
  return s[0];
}
static int SMAX=24;
static int decode(i128 N,int verbose){
  N0=N; found=0; nodes=0; totnodes=0;
  for(int S=0; S<=SMAX && !found; S++) chains(1,0,0,0,S);
  if(!found) return 0;
  i128 chk=verify(foundU,foundr);
  if(chk!=N){ char b1[64],b2[64]; p128(N,b1); p128(chk,b2); printf("VERIFY FAIL N=%s got %s\n",b1,b2); exit(1); }
  if(verbose){ char b[64]; p128(N,b); int k=0; for(i128 t=foundU;t;t>>=1) k++;
    printf("N=%s: source x=",b); for(int t=k-1;t>=0;t--) putchar(((foundU>>t)&1)?'a':'b'); if(k==0) printf("(empty)");
    printf("  readers:"); for(int q=foundr-1;q>=0;q--) printf(" c a^%d b^%d",seqi[q],seqj[q]);
    int nv=1+k; for(int q=0;q<foundr;q++) nv+=1+seqi[q]+seqj[q];
    printf("   [%d vertices, log2 N = %.1f, %llu solver nodes]\n",nv,(double)__builtin_log2((double)N),totnodes); }
  return 1;
}
static u64 rng=88172645463325252ULL; static u64 xr(void){ rng^=rng<<13; rng^=rng>>7; rng^=rng<<17; return rng; }
int main(int argc,char**argv){
  NODECAP=2000000;
  if(argc>1 && !strcmp(argv[1],"-r")){
    i128 LO=parse128(argv[2]), HI=parse128(argv[3]); int cnt=atoi(argv[4]); MAXR=argc>5?atoi(argv[5]):3; MAXIJ=argc>6?atoi(argv[6]):6;
    for(int i=0;i<=MAXIJ;i++) for(int j=0;j<=MAXIJ;j++) RD[i][j]=reader(i,j);
    int ok=0; clock_t c0=clock(); u64 totn=0;
    for(int t=0;t<cnt;t++){ i128 span=HI-LO; i128 r=((i128)xr()<<64)|xr(); if(r<0) r=-r; i128 N=LO+r%span; N|=1;
      int f=decode(N,0); totn+=totnodes; if(f) ok++; else { char b[64]; p128(N,b); printf("  no decode within caps: %s\n",b); } }
    printf("decoded %d of %d random odd targets in [LO,HI); %.3f s total, mean solver nodes %.0f\n",ok,cnt,(double)(clock()-c0)/CLOCKS_PER_SEC,(double)totn/cnt);
    return 0;
  }
  i128 N=parse128(argv[1]); MAXR=argc>2?atoi(argv[2]):3; MAXIJ=argc>3?atoi(argv[3]):6;
  for(int i=0;i<=MAXIJ;i++) for(int j=0;j<=MAXIJ;j++) RD[i][j]=reader(i,j);
  if(!decode(N,1)) printf("no decode within caps\n");
  return 0;
}
