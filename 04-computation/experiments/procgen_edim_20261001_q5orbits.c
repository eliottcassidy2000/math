// procgen_edim_20261001_q5orbits.c -- orbit representatives of a-subsets of V(Q5) under Aut(Q5) (order 3840).
// Canonical form = minimum 32-bit mask over the group; level a+1 = canonical forms of all one-point extensions
// of level-a reps (every (a+1)-set contains an a-subset, so this reaches every orbit). Counts are checked against
// Burnside and minimality is re-checked independently in procgen_edim_20261001_lib.check_q5_reps.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
static int perm[3840][32]; static int G=0;
static uint32_t T[3840][4][256];
static uint32_t img(int g,uint32_t m){return T[g][0][m&255]|T[g][1][(m>>8)&255]|T[g][2][(m>>16)&255]|T[g][3][m>>24];}
static uint32_t canon(uint32_t m){uint32_t best=0xffffffffu; for(int g=0;g<G;g++){uint32_t x=img(g,m); if(x<best)best=x;} return best;}
// hash set
static uint32_t *hs; static size_t hcap, hn;
static int hins(uint32_t k){ size_t h=((uint64_t)k*0x9E3779B97F4A7C15ULL)>>20; h&=hcap-1; while(hs[h]!=0xffffffffu){ if(hs[h]==k) return 0; h=(h+1)&(hcap-1);} hs[h]=k; hn++; return 1;}
int cmpu(const void*a,const void*b){uint32_t x=*(uint32_t*)a,y=*(uint32_t*)b;return x<y?-1:x>y;}
int main(int argc,char**argv){ int amax=atoi(argv[1]); const char*pre=argv[2];
  int p[5]={0,1,2,3,4}; // all 120 perms via lexicographic iteration
  int perms[120][5],np=0;
  for(int a=0;a<5;a++)for(int b=0;b<5;b++)for(int c=0;c<5;c++)for(int d=0;d<5;d++)for(int e=0;e<5;e++){
    int q[5]={a,b,c,d,e}; int used=0,ok=1; for(int i=0;i<5;i++){if(used>>q[i]&1)ok=0;used|=1<<q[i];} if(ok){memcpy(perms[np++],q,sizeof q);} }
  (void)p;
  for(int pi=0;pi<120;pi++) for(int t=0;t<32;t++){ for(int v=0;v<32;v++){ int w=v^t,u=0; for(int i=0;i<5;i++) if(w>>i&1) u|=1<<perms[pi][i]; perm[G][v]=u;} G++; }
  for(int g=0;g<G;g++) for(int bpos=0;bpos<4;bpos++) for(int val=0;val<256;val++){ uint32_t r=0; for(int j=0;j<8;j++) if(val>>j&1) r|=1u<<perm[g][bpos*8+j]; T[g][bpos][val]=r; }
  hcap=1<<22; hs=malloc(4*hcap);
  uint32_t *cur=malloc(4),*nxt; size_t ncur=1; cur[0]=0;
  for(int a=0;a<=amax;a++){
    char fn[256]; sprintf(fn,"%s_%d.txt",pre,a); FILE*f=fopen(fn,"w");
    for(size_t i=0;i<ncur;i++) fprintf(f,"%08x\n",cur[i]); fclose(f);
    printf("a=%d reps=%zu\n",a,ncur); fflush(stdout);
    if(a==amax) break;
    memset(hs,0xff,4*hcap); hn=0;
    for(size_t i=0;i<ncur;i++) for(int v=0;v<32;v++) if(!(cur[i]>>v&1)) hins(canon(cur[i]|(1u<<v)));
    nxt=malloc(4*hn); size_t j=0; for(size_t h=0;h<hcap;h++) if(hs[h]!=0xffffffffu) nxt[j++]=hs[h];
    qsort(nxt,j,4,cmpu); free(cur); cur=nxt; ncur=j;
  }
  return 0; }
