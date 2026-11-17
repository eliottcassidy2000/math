/* Wave 58: the triangle switch of three_visit_triangle_switch_20261117 with a transitive quotient and free shared-contact
   tapes on both sides; exact 3-segment counter (checked against Held-Karp), scans for (2^p+1)^2 and the value set. */
/* Triangle switch counter.  Exterior: transitive order of n vertices with the triangle module {b,c,y} (b->c->y->b)
   inserted at position t (after the first t exterior vertices).  Ports x,a: x->b, a->b, c->x, a->c, x->y, y->a.
   Exterior vertex i with bit g_i: g=1 beats both ports, g=0 beaten by both.  Minus: x->a.  Plus: a->x.
   H = sum over port orders (h1,h2) of 3-segment assignments S0 h1 S1 h2 S2 along the spine.
   Mode "check": compare with brute-force Held-Karp on random small cases.  Mode "scan N": all words, all positions. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
/* vertex ids for contacts: exterior bit g -> beats_x=g, beats_a=g.  triangle: b: beats_x=0,beats_a=0; c: beats_x=1, beats_a=0; y: beats_x=0, beats_a=1 */
enum {TB=0,TC=1,TY=2};
static const int tri_bx[3]={0,1,0}, tri_ba[3]={0,0,1};
static const int tri_arc[3][3]={{0,1,0},{0,0,1},{1,0,0}}; /* b->c, c->y, y->b */
/* port p in {0=x,1=a}; vertex beats port p? */
static int vb(int isTri,int id,int g,int p){ if(isTri) return p==0?tri_bx[id]:tri_ba[id]; return g; }
static u64 count_order(int n,int t,const int*g,int h1,int h2,int hub12){ /* hub12: 1 if h1->h2 */
  u64 dp[27]={0},nd[27]; dp[0]=1; /* st = s0 + 3 s1 + 9 s2 ; s: 0 empty, 1 ok, 2 bad */
  for(int pos=0;pos<=n;pos++){
    if(pos==t){ /* triangle macro */
      memset(nd,0,sizeof nd);
      for(int st=0;st<27;st++){ u64 c=dp[st]; if(!c) continue;
        for(int asg=0;asg<27;asg++){ int seg[3]={asg%3,asg/3%3,asg/9}; /* segment of b,c,y */
          /* for each segment collect its triangle vertices; enumerate path orders */
          int lst[3][3],cnt[3]={0,0,0}; for(int v=0;v<3;v++) lst[seg[v]][cnt[seg[v]]++]=v;
          /* orders per segment: list of (first,last) options with multiplicity */
          int nopt[3]; int fo[3][3],lo[3][3];
          int okall=1;
          for(int s=0;s<3;s++){ nopt[s]=0; if(cnt[s]==0){ nopt[s]=1; fo[s][0]=lo[s][0]=-1; continue; }
            if(cnt[s]==1){ fo[s][0]=lo[s][0]=lst[s][0]; nopt[s]=1; }
            else if(cnt[s]==2){ int u=lst[s][0],w=lst[s][1]; if(tri_arc[u][w]){fo[s][0]=u;lo[s][0]=w;} else {fo[s][0]=w;lo[s][0]=u;} nopt[s]=1; }
            else { /* rotations b c y, c y b, y b c */ int r[3][3]={{0,1,2},{1,2,0},{2,0,1}}; for(int k=0;k<3;k++){ fo[s][k]=r[k][0]; lo[s][k]=r[k][2]; } nopt[s]=3; } }
          if(!okall) continue;
          for(int o0=0;o0<nopt[0];o0++) for(int o1=0;o1<nopt[1];o1++) for(int o2=0;o2<nopt[2];o2++){
            int s0=st%3,s1=st/3%3,s2=st/9; int o[3]={o0,o1,o2}; int ok=1; int ns[3]={s0,s1,s2};
            for(int s=0;s<3&&ok;s++){ int f=fo[s][o[s]], l=lo[s][o[s]]; if(f<0) continue;
              int cur=ns[s];
              if(cur==0){ /* first vertex checks */
                if(s==1 && vb(1,f,0,h1)) ok=0;   /* need h1 -> f : f must NOT beat h1 */
                if(s==2 && vb(1,f,0,h2)) ok=0; }
              /* new last = l; ok status for right boundary */
              if(s==0) ns[s]= vb(1,l,0,h1)?1:2;   /* need l -> h1 */
              else if(s==1) ns[s]= vb(1,l,0,h2)?1:2;
              else ns[s]=1; }
            if(ok) nd[ns[0]+3*ns[1]+9*ns[2]]+=c; } } }
      memcpy(dp,nd,sizeof dp); }
    if(pos==n) break;
    int gg=g[pos]; memset(nd,0,sizeof nd);
    for(int st=0;st<27;st++){ u64 c=dp[st]; if(!c) continue; int s0=st%3,s1=st/3%3,s2=st/9;
      /* to S0 */ { int t0= gg?1:2; nd[t0+3*s1+9*s2]+=c; }
      /* to S1 */ if(s1!=0 || !gg){ int t1= gg?1:2; nd[s0+3*t1+9*s2]+=c; }
      /* to S2 */ if(s2!=0 || !gg){ nd[s0+3*s1+9*1]+=c; } }
    memcpy(dp,nd,sizeof dp); }
  u64 tot=0; for(int st=0;st<27;st++){ int s0=st%3,s1=st/3%3,s2=st/9; if(s0==2||s1==2||s2==2) continue; if(s1==0&&!hub12) continue; tot+=dp[st]; }
  return tot; }
static u64 Hsw(int n,int t,const int*g,int plus){ int xa = plus?0:1; /* x->a in minus */
  return count_order(n,t,g,0,1,xa)+count_order(n,t,g,1,0,!xa); }
/* brute force */
static u64 heldkarp(int N,const unsigned*out){ static u64 dp[1<<20][20]; for(unsigned S=0;S<(1u<<N);S++) for(int v=0;v<N;v++) dp[S][v]=0;
  for(int v=0;v<N;v++) dp[1u<<v][v]=1;
  for(unsigned S=1;S<(1u<<N);S++) for(int v=0;v<N;v++) if((S>>v&1)&&dp[S][v]) for(int w=0;w<N;w++) if(!(S>>w&1)&&(out[v]>>w&1)) dp[S|1u<<w][w]+=dp[S][v];
  u64 h=0; for(int v=0;v<N;v++) h+=dp[(1u<<N)-1][v]; return h; }
static void build(int n,int t,const int*g,int plus,unsigned*out,int*N){ /* vertices: exterior 0..n-1 in order but triangle inserted; ids: ext i -> i (spine order), b,c,y -> n,n+1,n+2, x->n+3, a->n+4 */
  int Nn=n+5; for(int i=0;i<Nn;i++) out[i]=0; int B=n,C=n+1,Y=n+2,X=n+3,A=n+4;
  for(int i=0;i<n;i++) for(int j=i+1;j<n;j++) out[i]|=1u<<j;
  for(int i=0;i<n;i++) for(int k=0;k<3;k++){ if(i<t) out[i]|=1u<<(n+k); else out[n+k]|=1u<<i; }
  out[B]|=1u<<C; out[C]|=1u<<Y; out[Y]|=1u<<B;
  out[X]|=1u<<B; out[A]|=1u<<B; out[C]|=1u<<X; out[A]|=1u<<C; out[X]|=1u<<Y; out[Y]|=1u<<A;
  for(int i=0;i<n;i++){ if(g[i]){ out[i]|=1u<<X; out[i]|=1u<<A; } else { out[X]|=1u<<i; out[A]|=1u<<i; } }
  if(plus) out[A]|=1u<<X; else out[X]|=1u<<A; *N=Nn; }

/* Modes: check | scan NMAX | values NMAX LIM */
int main(int argc,char**argv){
  if(argc>1 && !strcmp(argv[1],"check")){ srand(7); int bad=0,tot=0,gapbad=0; for(int it=0;it<400;it++){ int n=rand()%10, t=rand()%(n+1), g[32]; for(int i=0;i<n;i++) g[i]=rand()&1;
      for(int plus=0;plus<2;plus++){ unsigned out[32]; int N; build(n,t,g,plus,out,&N); u64 h1=heldkarp(N,out), h2=Hsw(n,t,g,plus); tot++; if(h1!=h2) bad++; }
      { int N; unsigned o1[32],o2[32]; build(n,t,g,0,o1,&N); build(n,t,g,1,o2,&N); if(heldkarp(N,o2)-heldkarp(N,o1)!=2) gapbad++; } }
    printf("triangle-switch counter vs Held-Karp: %d cases, mismatches %d; gap != 2 in %d of 400\n",tot,bad,gapbad); return 0; }
  if(argc>2 && !strcmp(argv[1],"scan")){ int NMAX=atoi(argv[2]); u64 T[20]; for(int p=1;p<20;p++){ u64 q=(1ULL<<p)+1; T[p]=q*q; } long hits=0, cfg=0;
    for(int n=0;n<=NMAX;n++) for(int t=0;t<=n;t++) for(u64 w=0;w<(1ULL<<n);w++){ int g[64]; for(int i=0;i<n;i++) g[i]=w>>i&1; u64 h=Hsw(n,t,g,0); cfg++;
        for(int p=1;p<20;p++) if(h==T[p]||h+2==T[p]){ hits++; printf("hit p=%d n=%d t=%d\n",p,n,t); } }
    printf("exterior <= %d: %ld configurations; minus or plus equal to (2^p+1)^2: %ld\n",NMAX,cfg,hits); return 0; }
  if(argc>3 && !strcmp(argv[1],"values")){ int NMAX=atoi(argv[2]); int LIM=atoi(argv[3]); unsigned char*hit=calloc(LIM+3,1);
    for(int n=0;n<=NMAX;n++) for(int t=0;t<=n;t++) for(u64 w=0;w<(1ULL<<n);w++){ int g[64]; for(int i=0;i<n;i++) g[i]=w>>i&1; u64 h=Hsw(n,t,g,0); if(h<=(u64)LIM) hit[h]=1; }
    int c=0; printf("minus values <= 200:"); for(int h=1;h<=200;h++) if(hit[h]) printf(" %d",h); printf("\n");
    for(int h=1;h<=LIM;h+=2) if(hit[h]) c++; printf("exterior <= %d: distinct odd minus values <= %d: %d of %d\n",NMAX,LIM,c,(LIM+1)/2); return 0; }
  return 0; }
