/* For every class rep of order <= 7 and every one-vertex extension (order 8) compute
   h = H(T), f = F_2(T) = sum_S H(T[S]) H(T[V\S]) (ordered two-path covers), strongness.
   Output unique (n,h,f,strong) lines; also first-appearance data for counts new at order 8. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static u64 dp[256][8], Hs[256];
static int nn; static int outm[8];
static void allsub(void){ int N=1<<nn; for(int S=0;S<N;S++){ for(int v=0;v<nn;v++) dp[S][v]=0; }
  for(int v=0;v<nn;v++) dp[1<<v][v]=1;
  for(int S=1;S<N;S++){ u64 tot=0; for(int v=0;v<nn;v++) if(S>>v&1){ if(S!=(1<<v)){ int P=S&~(1<<v); u64 s=0;
        for(int u=0;u<nn;u++) if((P>>u&1) && (outm[u]>>v&1)) s+=dp[P][u]; dp[S][v]=s; } tot+=dp[S][v]; } Hs[S]=tot; }
  Hs[0]=1; }
static int strongc(void){ int full=(1<<nn)-1; int r=1, st[8], sp=0; st[sp++]=0;
  while(sp){ int u=st[--sp]; int t=outm[u]&~r; while(t){ int w=__builtin_ctz(t); t&=t-1; r|=1<<w; st[sp++]=w; } }
  if(r!=full) return 0; r=1; sp=0; st[sp++]=0;
  while(sp){ int u=st[--sp]; int in=0; for(int w=0;w<nn;w++) if(outm[w]>>u&1) in|=1<<w; int t=in&~r; while(t){ int w=__builtin_ctz(t); t&=t-1; r|=1<<w; st[sp++]=w; } }
  return r==full; }
typedef struct { int n; u64 h, f; int s; } rec;
static rec *R; static size_t nr=0, cap=0;
static void push(int n,u64 h,u64 f,int s){ if(nr==cap){ cap=cap?cap*2:1024; R=realloc(R,cap*sizeof(rec)); } R[nr++]=(rec){n,h,f,s}; }
static int cmp(const void*a,const void*b){ const rec*x=a,*y=b; if(x->n!=y->n) return x->n-y->n; if(x->h!=y->h) return x->h<y->h?-1:1; if(x->f!=y->f) return x->f<y->f?-1:1; return x->s-y->s; }
static char seenH[1<<20];
int main(void){ FILE*fi=fopen("w57_reps.txt","r"); int reps7[456][7], k7=0; char line[512];
  while(fgets(line,sizeof line,fi)){ char*p=line; int n=strtol(p,&p,10); nn=n; for(int i=0;i<n;i++) outm[i]=strtol(p,&p,10);
    allsub(); u64 h=Hs[(1<<n)-1], f=0; for(int S=0;S<(1<<n);S++) f+=Hs[S]*Hs[((1<<n)-1)&~S];
    int s=strongc(); push(n,h,f,s); if(h<(1u<<20)) seenH[h]=1; if(n==7){ memcpy(reps7[k7++],outm,sizeof(int)*7); } }
  fclose(fi);
  /* order 8 */
  static u64 parmask[1<<20][4]; static unsigned char newstrong[1<<20]; /* bit0 some strong, bit1 some nonstrong */
  for(int r=0;r<k7;r++) for(int S=0;S<128;S++){ nn=8; for(int i=0;i<7;i++){ outm[i]=reps7[r][i]; if(!(S>>i&1)) outm[i]|=1<<7; } outm[7]=S;
      allsub(); u64 h=Hs[255], f=0; for(int T=0;T<256;T++) f+=Hs[T]*Hs[255&~T]; int s=strongc(); push(8,h,f,s);
      if(h<(1u<<20) && !seenH[h]){ newstrong[h]|= s?1:2; for(int v=0;v<8;v++){ u64 p=Hs[255&~(1<<v)]; if(p<256) parmask[h][p>>6]|=1ULL<<(p&63); else parmask[h][3]|=1ULL<<63; } } }
  qsort(R,nr,sizeof(rec),cmp); FILE*fo=fopen("w57_pairs.txt","w"); size_t u=0;
  for(size_t i=0;i<nr;i++){ if(i&&!cmp(&R[i],&R[i-1])) continue; fprintf(fo,"%d %llu %llu %d\n",R[i].n,R[i].h,R[i].f,R[i].s); u++; } fclose(fo);
  fprintf(stderr,"unique (n,h,f,strong) records: %zu\n",u);
  FILE*fn=fopen("w57_first8.txt","w"); for(int h=1;h<(1<<20);h++) if(newstrong[h]){ fprintf(fn,"%d %d",h,newstrong[h]);
      for(int p=0;p<255;p++) if(parmask[h][p>>6]>>(p&63)&1) fprintf(fn," %d",p); if(parmask[h][3]>>63&1) fprintf(fn," big"); fprintf(fn,"\n"); }
  fclose(fn); return 0; }
