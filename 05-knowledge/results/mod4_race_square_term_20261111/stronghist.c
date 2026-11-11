// stronghist.c - histogram of H over the STRONG tournament classes of order n (and over all classes).
// stdin: gentourng output (upper triangle row by row, '1' = i->j for i<j).  Usage: stronghist n
// strong  <=>  for the scores sorted increasingly, s_1+...+s_k > C(k,2) for every 1 <= k <= n-1 (Moon / Landau).
// Output: one line "h strong all" per odd h that occurs, then a summary line on stderr.
#include <stdio.h>
#include <stdlib.h>
typedef unsigned long long u64;
#define MAXN 11
#define HMAX 200000
static u64 f[1<<MAXN][MAXN]; static unsigned hs[HMAX], ha[HMAX];
int main(int argc,char**argv){ int n=atoi(argv[1]); int full=(1<<n)-1; char line[256]; int inm[MAXN], sc[MAXN];
  u64 nt=0, ns=0, over=0;
  while(fgets(line,sizeof line,stdin)){ if(line[0]!='0'&&line[0]!='1')continue; int k=0;
    for(int i=0;i<n;i++){inm[i]=0;sc[i]=0;}
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){ if(line[k++]=='1'){inm[j]|=1<<i;sc[i]++;} else {inm[i]|=1<<j;sc[j]++;} }
    nt++;
    for(int S=1;S<=full;S++) for(int v=0;v<n;v++){ if(!(S>>v&1)){f[S][v]=0;continue;} int R=S&~(1<<v);
        if(R==0){f[S][v]=1;continue;} u64 a=0; int mm=R&inm[v]; while(mm){int w=__builtin_ctz(mm);mm&=mm-1;a+=f[R][w];} f[S][v]=a; }
    u64 H=0; for(int v=0;v<n;v++) H+=f[full][v];
    // strong test
    int s[MAXN]; for(int i=0;i<n;i++) s[i]=sc[i];
    for(int i=1;i<n;i++){ int t=s[i], j=i-1; while(j>=0&&s[j]>t){s[j+1]=s[j];j--;} s[j+1]=t; }
    int strong=1, acc=0; for(int kk=1;kk<n;kk++){ acc+=s[kk-1]; if(acc<=kk*(kk-1)/2){strong=0;break;} }
    if(n==1) strong=1;
    if(H>=HMAX){over++;continue;}
    ha[H]++; if(strong){hs[H]++;ns++;}
  }
  for(int h=1;h<HMAX;h+=2) if(ha[h]) printf("%d %u %u\n", h, hs[h], ha[h]);
  fprintf(stderr,"order %d: classes %llu, strong %llu, overflow %llu\n", n, nt, ns, over);
  return 0; }
