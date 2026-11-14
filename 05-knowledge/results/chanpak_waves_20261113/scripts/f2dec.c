// f2dec: target-first decoder for root-change (O/Z/J) words with at most RMAX root changes and ARBITRARY reader
// blocks J w (w any O/Z word): the family F_{<=2} of HYP-9494 (D) when RMAX = 2.
// Letters: a = O, b = Z, c = J.  Root triple recurrence as in catdec.c.  The readers are enumerated by total length,
// the source (one-apex word, any length) is solved low bit first with the readers' factors carried inside
// (A g(U) + B U = T), and every output word is re-verified by the 4D recurrence.
// usage: f2dec -r LO HI count RMAX MAXLEN     (random odd targets)    |   f2dec N RMAX MAXLEN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
typedef __int128 i128; typedef unsigned long long u64;
static i128 Hc(i128 U){ if(U==0) return 1; int k=0; for(i128 t=U;t;t>>=1) k++; i128 s=1;
  for(int t=0;t<k;t++){ int bt=(int)((U>>t)&1); if(!bt) s+=U>>(t+1); if(t+1<k && ((U>>(t+1))&1) && !bt) s+=1; } return s; }
static void p128(i128 x,char*buf){ char t[64]; int n=0; if(x==0) t[n++]='0'; while(x){t[n++]='0'+(int)(x%10); x/=10;} int m=0; while(n) buf[m++]=t[--n]; buf[m]=0; }
static i128 parse128(const char*s){ i128 x=0; while(*s){ x=x*10+(*s-'0'); s++;} return x; }
static u64 nodes, NODECAP, totnodes;
static i128 solveU(i128 A,i128 B,i128 T,int req,int depth){
  nodes++; if(nodes>NODECAP) return -1;
  if(T<A) return -1;
  if(T==A && req!=2) return 0;
  if(B>0){ i128 lo=(T-A)/(A+B); if(lo<1) lo=1; i128 hi=(T-A)/B; if(lo>hi) return -1;
    if(hi-lo<=8){ for(i128 U=lo;U<=hi;U++){ if(req==1&&(U&1)) continue; if(req==2&&!(U&1)) continue; if(A*Hc(U)+B*U==T) return U; } return -1; } }
  if(depth>125) return -1;
  if(req!=1 && B>0){ i128 r=solveU(A,2*B,T-B,0,depth+1); if(r>=0) return 2*r+1; }
  if(req!=2){ i128 r=solveU(A,2*B+A,T-A,2,depth+1); if(r>=0) return 2*r; r=solveU(A,2*B+A,T,1,depth+1); if(r>=1) return 2*r; }
  return -1;
}
static void st_a(i128*s){ s[1]=0; s[2]=2*s[2]; }
static void st_b(i128*s){ i128 H=s[0],E=s[1],L=s[2],h=s[3]; s[0]=H+L-E; s[1]=h; s[2]=2*L-h; }
static void st_c(i128*s){ i128 H=s[0],E=s[1],L=s[2]; s[0]=H+L-E; s[1]=H-E; s[2]=H+L; s[3]=H; }
// block J w: w given by (bits,len), bit 1 = O (a), 0 = Z (b), first letter = most significant bit
static void blockmat(u64 bits,int len,i128 M[3][3]){
  for(int col=0;col<3;col++){ i128 s[4]={0,0,0,0}; s[col]=1; st_c(s);
    for(int t=len-1;t>=0;t--){ if((bits>>t)&1) st_a(s); else st_b(s); }
    for(int row=0;row<3;row++) M[row][col]=s[row]; }
}
static i128 N0; static int RMAX, MAXLEN; static int found; static i128 fU; static u64 fb[4]; static int fl[4], fr;
static u64 cb[4]; static int cl[4];
static int try_source(i128 A,i128 C,i128 B){
  if(A<=0||B<0) return 0; if(A+B+C>N0) return 0;
  double est=sqrt((double)N0*(double)A)/(double)(B>0?B:1);
  nodes=0; NODECAP=(u64)(6*est)+512; if(NODECAP>200000000ULL) NODECAP=200000000ULL;
  i128 U=solveU(A,B,N0-B,2,0); if(U<0) U=solveU(A,B,N0-B-C,1,0);
  totnodes+=nodes; if(U>=0){ fU=U; return 1; } return 0;
}
// enumerate reader chains with total letter count exactly S (r readers), peeled from the end; phi = e1 * M_last * ...
static void chains(i128 A,i128 C,i128 B,int r,int S){
  if(found) return; if(A+B+C>N0) return;
  if(S==0){ if(try_source(A,C,B)){ found=1; fr=r; memcpy(fb,cb,sizeof fb); memcpy(fl,cl,sizeof fl);} }
  if(r>=RMAX) return;
  for(int len=0; len<=S && len<=MAXLEN && !found; len++){
    if(S==0 && len>0) break;
    for(u64 bits=0; bits<(1ULL<<len) && !found; bits++){
      i128 M[3][3]; blockmat(bits,len,M);
      i128 nA=A*M[0][0]+C*M[1][0]+B*M[2][0], nC=A*M[0][1]+C*M[1][1]+B*M[2][1], nB=A*M[0][2]+C*M[1][2]+B*M[2][2];
      cb[r]=bits; cl[r]=len;
      chains(nA,nC,nB,r+1,S-len);
    }
  }
}
static i128 verify(void){
  i128 s[4]={1,1,1,1}; int k=0; for(i128 t=fU;t;t>>=1) k++;
  for(int t=k-1;t>=0;t--){ if((fU>>t)&1) st_a(s); else st_b(s); }
  for(int q=fr-1;q>=0;q--){ st_c(s); for(int t=fl[q]-1;t>=0;t--){ if((fb[q]>>t)&1) st_a(s); else st_b(s);} }
  return s[0];
}
static int decode(i128 N,int verbose){
  N0=N; found=0; totnodes=0;
  for(int S=0; S<=RMAX*MAXLEN && !found; S++) chains(1,0,0,0,S);
  if(!found) return 0;
  if(verify()!=N){ char b[64]; p128(N,b); printf("VERIFY FAIL %s\n",b); exit(1); }
  if(verbose){ char b[64]; p128(N,b); int k=0; for(i128 t=fU;t;t>>=1) k++;
    printf("N=%s: source ",b); for(int t=k-1;t>=0;t--) putchar(((fU>>t)&1)?'O':'Z'); if(!k) printf("(empty)");
    for(int q=fr-1;q>=0;q--){ printf(" J"); for(int t=fl[q]-1;t>=0;t--) putchar(((fb[q]>>t)&1)?'O':'Z'); }
    printf("   [%d root changes, %llu solver nodes]\n",fr,totnodes); }
  return 1;
}
static u64 rng=0x9E3779B97F4A7C15ULL; static u64 xr(void){ rng^=rng<<13; rng^=rng>>7; rng^=rng<<17; return rng; }
int main(int argc,char**argv){
  if(!strcmp(argv[1],"-r")){
    i128 LO=parse128(argv[2]),HI=parse128(argv[3]); int cnt=atoi(argv[4]); RMAX=atoi(argv[5]); MAXLEN=atoi(argv[6]);
    if(argc>7) rng^=strtoull(argv[7],0,10)*0x2545F4914F6CDD1DULL;
    int ok=0; clock_t c0=clock(); double tn=0;
    for(int t=0;t<cnt;t++){ i128 r=((i128)xr()<<64)|xr(); if(r<0) r=-r; i128 N=LO+r%(HI-LO); N|=1;
      if(decode(N,1)) ok++; else { char b[64]; p128(N,b); printf("  NOT DECODED (caps): %s\n",b); } tn+=totnodes; fflush(stdout); }
    printf("decoded %d of %d; %.2f s total; mean nodes %.3g\n",ok,cnt,(double)(clock()-c0)/CLOCKS_PER_SEC,tn/cnt); return 0; }
  i128 N=parse128(argv[1]); RMAX=atoi(argv[2]); MAXLEN=atoi(argv[3]);
  if(!decode(N,1)) printf("not decoded within caps\n"); return 0;
}
