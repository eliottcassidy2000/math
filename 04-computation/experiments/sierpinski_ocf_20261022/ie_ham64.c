/* ie_ham64.c (64-bit vertex masks, circulant mode for n <= 63; copy of paley_extremal_20261021/ie_ham.c) -- exact Hamiltonian-path counts by inclusion-exclusion over vertex subsets,
 *   H = sum_{S subseteq V} (-1)^{n-|S|} W(S),  W(S) = #walks with n vertices inside S,
 * computed in Z/2^128 (exact since H < 2^128 for n <= 33).  Symmetry reduction:
 *   circ mode: circulant Cay(Z_n, conn); subsets up to rotation = binary necklaces (FKM), weight = period.
 *   f27 mode : Paley P_27 on F_27 = F_3[t]/(t^3 - t - 1); subsets up to translation by F_27 (canonical min form).
 * mac-mini-2026-10-21-fanopsp.  Usage: ./ie_ham circ n s1 s2 ...   |   ./ie_ham f27 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned __int128 u128; typedef uint64_t u64;
static int n; static uint64_t out[64];
static u128 W(uint64_t S){ u128 v[40], w[40]; int idx[40], k=0; for(int i=0;i<n;i++) if(S>>i&1){ idx[k++]=i; }
  if(k==0) return 0; for(int a=0;a<k;a++) v[idx[a]]=1;
  for(int step=1;step<n;step++){ for(int a=0;a<k;a++){ int i=idx[a]; uint64_t m=out[i]&S; u128 s=0; while(m){ int j=__builtin_ctzll(m); m&=m-1; s+=v[j]; } w[i]=s; }
    for(int a=0;a<k;a++) v[idx[a]]=w[idx[a]]; }
  u128 tot=0; for(int a=0;a<k;a++) tot+=v[idx[a]]; return tot; }
static void p128(u128 x){ char b[64]; int i=63; b[i]=0; if(!x){printf("0");return;} while(x){ b[--i]='0'+(int)(x%10); x/=10; } printf("%s",b+i); }
/* --- necklaces --- */
static int a[64]; static u128 acc; static long nneck;
static void visit(int p){ uint64_t S=0; int sz=0; for(int t=1;t<=n;t++) if(a[t]){ S|=1ull<<(t-1); sz++; }
  u128 w=W(S)*(u128)p; if((n-sz)&1) acc-=w; else acc+=w; nneck++; }
static void gen(int t,int p){ if(t>n){ if(n%p==0) visit(p); return; } a[t]=a[t-p]; gen(t+1,p); if(a[t-p]==0){ a[t]=1; gen(t+1,t); } }
int main(int argc,char**argv){
  if(!strcmp(argv[1],"circ")){ n=atoi(argv[2]); int inS[64]={0}; for(int i=3;i<argc;i++) inS[atoi(argv[i])%n]=1;
    for(int i=0;i<n;i++){ out[i]=0; for(int j=0;j<n;j++) if(j!=i && inS[((j-i)%n+n)%n]) out[i]|=1ull<<j; }
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++) if(((out[i]>>j)&1ull)==((out[j]>>i)&1ull)){ fprintf(stderr,"not a tournament\n"); return 1; }
    acc=0; nneck=0; a[0]=0; gen(1,1);
    printf("circ n=%d necklaces=%ld H=",n,nneck); p128(acc); printf("\n"); return 0; }
  if(!strcmp(argv[1],"f27")){ n=27; /* element = c0 + c1 t + c2 t^2, index c0+3c1+9c2; t^3 = t + 1 */
    int add[27][27], mul[27][27];
    for(int x=0;x<27;x++) for(int y=0;y<27;y++){ int xs[3]={x%3,(x/3)%3,x/9}, ys[3]={y%3,(y/3)%3,y/9}; int s[3]; for(int i=0;i<3;i++) s[i]=(xs[i]+ys[i])%3; add[x][y]=s[0]+3*s[1]+9*s[2];
      int pr[5]={0}; for(int i=0;i<3;i++) for(int j=0;j<3;j++) pr[i+j]=(pr[i+j]+xs[i]*ys[j])%3;
      /* reduce t^4 = t^2 + t, t^3 = t + 1 */
      pr[2]=(pr[2]+pr[4])%3; pr[1]=(pr[1]+pr[4])%3; pr[4]=0; pr[1]=(pr[1]+pr[3])%3; pr[0]=(pr[0]+pr[3])%3; pr[3]=0;
      mul[x][y]=pr[0]+3*pr[1]+9*pr[2]; }
    int sq[27]={0}; for(int x=1;x<27;x++) sq[mul[x][x]]=1; int nsq=0; for(int x=1;x<27;x++) nsq+=sq[x];
    int neg[27]; for(int x=0;x<27;x++){ int xs[3]={x%3,(x/3)%3,x/9}; neg[x]=((3-xs[0])%3)+3*((3-xs[1])%3)+9*((3-xs[2])%3); }
    for(int x=1;x<27;x++) if(sq[x]==sq[neg[x]]){ fprintf(stderr,"-1 is a square? not a tournament\n"); return 1; }
    for(int i=0;i<27;i++){ out[i]=0; for(int j=0;j<27;j++) if(j!=i && sq[add[j][neg[i]]]) out[i]|=1u<<j; }
    /* translation tables: image of each vertex under +b */
    static uint32_t img[27][4][256]; for(int b=0;b<27;b++) for(int byte=0;byte<4;byte++) for(int val=0;val<256;val++){ uint32_t r=0; for(int bit=0;bit<8;bit++) if(val>>bit&1){ int x=byte*8+bit; if(x<27) r|=1u<<add[x][b]; } img[b][byte][val]=r; }
    u128 tot=0; long reps=0; uint32_t FULL=(1u<<27)-1;
    for(uint32_t S=0;S<=FULL;S++){ /* canonical iff S <= all translates; weight = orbit size = 27/|stab| */
      int canon=1, stab=0; for(int b=0;b<27;b++){ uint32_t T=img[b][0][S&255]|img[b][1][(S>>8)&255]|img[b][2][(S>>16)&255]|img[b][3][(S>>24)&255]; if(T<S){canon=0;break;} if(T==S) stab++; }
      if(!canon) continue; reps++; int sz=__builtin_popcount(S); u128 w=W(S)*(u128)(27/stab); if((27-sz)&1) tot-=w; else tot+=w; }
    printf("P27 squares=%d orbit-reps=%ld H=",nsq,reps); p128(tot); printf("\n"); return 0; }
  return 1; }
