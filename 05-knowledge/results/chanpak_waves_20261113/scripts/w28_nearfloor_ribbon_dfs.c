/* nfrib.c -- near-floor span-w ribbon words by DFS with a level-dependent prune.
   Records, for each length k <= M, every strong word with H < Cf * f(k): value bitmaps, counts,
   min H coprime to 15, min prime H.  Prefix prune: H(prefix of length k) < Cp * f(k).
   usage: nfrib matrices automaton w M Cp Cf outprefix */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define MAXD 128
#define MAXL 8
#define MAXC 256
#define MAXM 40
static int d, nl, nf, w, M, nc, na;
static int fin[MAXD], lmask[MAXL], nnz[MAXL];
static int *Mj[MAXL], *Mi[MAXL];
static uint64_t *Mm[MAXL];
static int trans[MAXC][MAXL];
static uint8_t accept[MAXC];
static double Cp, Cf;
static uint64_t fl[MAXM+1], lim[MAXM+1], flim[MAXM+1];
static uint8_t *bits[MAXM+1];
static uint64_t nwords[MAXM+1], nodes=0, minc15[MAXM+1], minpr[MAXM+1];
static uint64_t wit15[MAXM+1][2], witpr[MAXM+1][2]; static uint64_t minall[MAXM+1], minc3[MAXM+1], minc5[MAXM+1];
static int isprime(uint64_t n){ if(n<2) return 0; for(uint64_t p=2;p*p<=n;p++) if(n%p==0) return 0; return 1; }
static uint64_t floorf_(int n){ uint64_t best=0; for(int b=0;3*b<=n-1;b++){int r=n-1-3*b; if(r%2) continue; int a=r/2; uint64_t v=1; for(int i=0;i<a;i++) v*=3; for(int i=0;i<b;i++) v*=5; if(!best||v<best) best=v;} return best; }
static inline uint64_t Hval(const uint64_t *v){ uint64_t s=0; for(int k=0;k<nf;k++) s+=v[fin[k]]; return s; }
static uint64_t pk[2];
static void dfs(const uint64_t *v, int depth, int cs){
  nodes++;
  if(depth>=M) return;
  int q=depth+1;
  for(int li=0;li<nl;li++){
    int m=lmask[li];
    if(q>=2 && (m>>(q-1))) continue;
    if(q==1 && m) continue;
    int c2=trans[cs][li]; if(c2<0) continue;
    uint64_t nv[MAXD]; memset(nv,0,sizeof(uint64_t)*d);
    const int *J=Mj[li],*I=Mi[li]; const uint64_t *MM=Mm[li];
    for(int k=0;k<nnz[li];k++) nv[J[k]]+=MM[k]*v[I[k]];
    uint64_t H=Hval(nv);
    if(H>=lim[q]) continue;
    uint64_t save0=pk[0],save1=pk[1];
    int bitpos=w*depth; if(bitpos<64) pk[0]|=((uint64_t)m)<<bitpos; else pk[1]|=((uint64_t)m)<<(bitpos-64);
    if(accept[c2] && H<flim[q]){
      nwords[q]++; if(!minall[q]||H<minall[q]) minall[q]=H; if(H%3 && (!minc3[q]||H<minc3[q])) minc3[q]=H; if(H%5 && (!minc5[q]||H<minc5[q])) minc5[q]=H;
      uint64_t idx=H>>1; bits[q][idx>>3]|=(uint8_t)(1u<<(idx&7));
      if(H%3 && H%5 && (!minc15[q]||H<minc15[q])){minc15[q]=H; wit15[q][0]=pk[0]; wit15[q][1]=pk[1];}
      if((!minpr[q]||H<minpr[q]) && isprime(H)){minpr[q]=H; witpr[q][0]=pk[0]; witpr[q][1]=pk[1];}
    }
    dfs(nv,depth+1,c2);
    pk[0]=save0; pk[1]=save1;
  }
}
static void pword(FILE *o, uint64_t a, uint64_t b, int len){
  for(int q=0;q<len;q++){ int bp=w*q; int m= bp<64 ? (int)((a>>bp)&((1u<<w)-1)) : (int)((b>>(bp-64))&((1u<<w)-1)); fprintf(o,"%s%d",q?".":"",m);} }
int main(int argc,char**argv){
  if(argc<8){fprintf(stderr,"usage\n");return 1;}
  FILE*f=fopen(argv[1],"r"); w=atoi(argv[3]); M=atoi(argv[4]); Cp=atof(argv[5]); Cf=atof(argv[6]);
  if(fscanf(f,"%d %d %d",&d,&nl,&nf)!=3) return 1;
  for(int k=0;k<nf;k++) if(fscanf(f,"%d",&fin[k])!=1) return 1;
  for(int li=0;li<nl;li++){ if(fscanf(f,"%d %d",&lmask[li],&nnz[li])!=2) return 1;
    Mj[li]=malloc(sizeof(int)*nnz[li]); Mi[li]=malloc(sizeof(int)*nnz[li]); Mm[li]=malloc(sizeof(uint64_t)*nnz[li]);
    for(int k=0;k<nnz[li];k++){unsigned long long mm; if(fscanf(f,"%d %d %llu",&Mj[li][k],&Mi[li][k],&mm)!=3) return 1; Mm[li][k]=mm;} }
  fclose(f);
  f=fopen(argv[2],"r"); if(fscanf(f,"%d %d",&nc,&na)!=2) return 1;
  for(int k=0;k<na;k++){int s; if(fscanf(f,"%d",&s)!=1) return 1; accept[s]=1;}
  for(int s=0;s<nc;s++) for(int li=0;li<nl;li++) if(fscanf(f,"%d",&trans[s][li])!=1) return 1;
  fclose(f);
  for(int k=1;k<=M;k++){ fl[k]=floorf_(k); if(k<3) fl[k]=1;} for(int k=1;k<=M;k++){ lim[k]=(uint64_t)(Cp*fl[k])+1; if(k<=6) lim[k]=(uint64_t)(Cp*fl[6]*3)+1; flim[k]=(uint64_t)(Cf*fl[k]); bits[k]=calloc(flim[k]/16+2,1);}
  uint64_t v0[MAXD]; memset(v0,0,sizeof v0); v0[0]=1;
  dfs(v0,0,0);
  char name[512]; snprintf(name,sizeof name,"%s.txt",argv[7]); FILE*o=fopen(name,"w");
  fprintf(o,"# span %d, M %d, Cp %.2f, Cf %.2f, nodes %llu\n",w,M,Cp,Cf,(unsigned long long)nodes);
  for(int k=3;k<=M;k++){
    uint64_t F=fl[k]; uint64_t att=0, slots=0, att2=0, slots2=0; uint64_t lo2=(5*F+2)/3; /* window [5f/3,3f) */
    for(uint64_t N=F|1; N<flim[k]; N+=2){ slots++; uint64_t idx=N>>1; int a=(bits[k][idx>>3]>>(idx&7))&1; att+=a; if(N>=lo2 && N<3*F){slots2++; att2+=a;} }
    fprintf(o,"k %d f %llu words %llu | [f,%.1ff) odd attained %llu/%llu | window [5f/3,3f) %llu/%llu | min coprime15 %llu (%.3f f) | min prime %llu (%.3f f)\n",
      k,(unsigned long long)F,(unsigned long long)nwords[k],Cf,(unsigned long long)att,(unsigned long long)slots,(unsigned long long)att2,(unsigned long long)slots2,
      (unsigned long long)minc15[k], minc15[k]?(double)minc15[k]/F:0,(unsigned long long)minpr[k], minpr[k]?(double)minpr[k]/F:0);
    fprintf(o,"   min %llu (%.3f) min not-div-3 %llu (%.3f) min not-div-5 %llu (%.3f)\n",(unsigned long long)minall[k],(double)minall[k]/F,(unsigned long long)minc3[k],(double)minc3[k]/F,(unsigned long long)minc5[k],(double)minc5[k]/F); fprintf(o,"   wit15 "); pword(o,wit15[k][0],wit15[k][1],k); fprintf(o,"   witpr "); pword(o,witpr[k][0],witpr[k][1],k); fprintf(o,"\n");
  }
  /* list misses in window for each k */
  for(int k=7;k<=M;k++){ uint64_t F=fl[k]; uint64_t lo2=(5*F+2)/3; fprintf(o,"misses k=%d in [5f/3,3f):",k); int cnt=0;
    for(uint64_t N=lo2|1; N<3*F; N+=2){ uint64_t idx=N>>1; if(!((bits[k][idx>>3]>>(idx&7))&1)){ if(cnt<40) fprintf(o," %llu",(unsigned long long)N); cnt++;} }
    fprintf(o," (total %d)\n",cnt); }
  fclose(o);
  printf("nodes %llu\n",(unsigned long long)nodes);
  return 0;
}
