/* Wave 59: all classes of orders 1..8 (order 8 by canonical forms of one-vertex extensions) -> w59_reps8.txt */
/* Order-8 classes by canonical form (extensions of the 456 order-7 reps), then all one-vertex extensions to order 9.
   Reports: #classes at 8 (expect 6880), spectrum at 9, new counts at 9 with strongness of realizers and deletion parents. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static u64 dp[512][9], Hs[512]; static int nn, outm[9];
static void allsub(void){ int N=1<<nn; for(int S=0;S<N;S++) for(int v=0;v<nn;v++) dp[S][v]=0;
  for(int v=0;v<nn;v++) dp[1<<v][v]=1;
  for(int S=1;S<N;S++){ u64 tot=0; for(int v=0;v<nn;v++) if(S>>v&1){ if(S!=(1<<v)){ int P=S&~(1<<v); u64 s=0;
        for(int u=0;u<nn;u++) if((P>>u&1)&&(outm[u]>>v&1)) s+=dp[P][u]; dp[S][v]=s; } tot+=dp[S][v]; } Hs[S]=tot; } Hs[0]=1; }
static int strongc(void){ int full=(1<<nn)-1, r=1, st[9], sp=0; st[sp++]=0;
  while(sp){ int u=st[--sp]; int t=outm[u]&~r; while(t){ int w=__builtin_ctz(t); t&=t-1; r|=1<<w; st[sp++]=w; } }
  if(r!=full) return 0; r=1; sp=0; st[sp++]=0;
  while(sp){ int u=st[--sp]; int in=0; for(int w=0;w<nn;w++) if(outm[w]>>u&1) in|=1<<w; int t=in&~r; while(t){ int w=__builtin_ctz(t); t&=t-1; r|=1<<w; st[sp++]=w; } }
  return r==full; }
/* canonical form for n=8: invariant = (score, H(T-v)); permute within cells */
static int perm[8], best_set; static unsigned bestcode; static int cellstart[9], ncell, order[8];
static unsigned code_of(const int*p){ /* p[i] = original vertex placed at position i */ unsigned c=0; int t=0;
  for(int i=0;i<8;i++) for(int j=i+1;j<8;j++){ if(outm[p[i]]>>p[j]&1) c|=1u<<t; t++; } return c; }
static void rec(int pos){ if(pos==8){ unsigned c=code_of(perm); if(!best_set||c<bestcode){ bestcode=c; best_set=1; } return; }
  /* find cell of pos */ int ce=0; while(!(cellstart[ce]<=pos && pos<cellstart[ce+1])) ce++;
  for(int k=cellstart[ce];k<cellstart[ce+1];k++){ int v=order[k]; int used=0; for(int q=cellstart[ce];q<pos;q++) if(perm[q]==v) used=1; if(used) continue; perm[pos]=v; rec(pos+1); } }
static long long inv[8];
static int cmpv(const void*a,const void*b){ long long x=inv[*(const int*)a], y=inv[*(const int*)b]; return x<y?-1:x>y; }
static unsigned canon8(void){ nn=8; allsub(); for(int v=0;v<8;v++){ int sc=__builtin_popcount(outm[v]); inv[v]=(long long)sc*100000+(long long)Hs[255&~(1<<v)]; order[v]=v; }
  qsort(order,8,sizeof(int),cmpv); ncell=0; cellstart[0]=0; for(int i=1;i<8;i++) if(inv[order[i]]!=inv[order[i-1]]) cellstart[++ncell]=i; cellstart[++ncell]=8;
  best_set=0; rec(0); return bestcode; }

/* dump all classes of orders 1..8: orders <= 7 from w57_reps.txt, order 8 by canonical forms of one-vertex extensions */
int main(void){ FILE*fi=fopen("w57_reps.txt","r"); FILE*fo=fopen("w59_reps8.txt","w"); int reps7[456][7], k7=0; char line[512];
  while(fgets(line,sizeof line,fi)){ fputs(line,fo); char*p=line; int n=strtol(p,&p,10); int o[8]; for(int i=0;i<n;i++) o[i]=strtol(p,&p,10); if(n==7) memcpy(reps7[k7++],o,sizeof(int)*7); }
  fclose(fi); unsigned char *seen=calloc(1u<<25,1); int k8=0;
  for(int r=0;r<k7;r++) for(int S=0;S<128;S++){ for(int i=0;i<7;i++){ outm[i]=reps7[r][i]; if(!(S>>i&1)) outm[i]|=1<<7; } outm[7]=S;
      unsigned c=canon8(); if(!(seen[c>>3]>>(c&7)&1)){ seen[c>>3]|=1<<(c&7); k8++; fprintf(fo,"8"); for(int i=0;i<8;i++) fprintf(fo," %d",outm[i]); fprintf(fo,"\n"); } }
  fclose(fo); fprintf(stderr,"order-8 classes: %d\n",k8); return 0; }
