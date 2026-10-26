/* For each class: enumerate the labelled C3-class (same cyclic triples) and count distinct score FUNCTIONS among its members,
   plus report whether the class size equals the number of distinct score functions (i.e. (C3, score function) determines T). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int N; static int target[12][12][12]; static int arc[12][12]; static int pairs[66][2], np;
static long long cnt; static unsigned long long scores[20000]; static int ns;
static int cyc(int a,int b,int c){ int ab=a<b?arc[a][b]:!arc[b][a]; int bc=b<c?arc[b][c]:!arc[c][b]; int ca=c<a?arc[c][a]:!arc[a][c]; return (ab&&bc&&ca)||(!ab&&!bc&&!ca); }
static void dfs(int p){
    if(p==np){ cnt++; unsigned long long key=0; for(int i=0;i<N;i++){int s=0; for(int j=0;j<N;j++) if(j!=i) s+= (i<j?arc[i][j]:!arc[j][i]); key=key*16+s;} if(ns<20000) scores[ns++]=key; return; }
    int i=pairs[p][0], j=pairs[p][1];
    for(int d=0;d<2;d++){ arc[i][j]=d; int ok=1;
        for(int k=0;k<N&&ok;k++){ if(k==i||k==j) continue; int a=i<k?i:k,b=i<k?k:i,c=j<k?j:k,e=j<k?k:j; if(arc[a][b]<0||arc[c][e]<0) continue;
            int lo=i,mi=j,hi=k; if(lo>mi){int t=lo;lo=mi;mi=t;} if(mi>hi){int t=mi;mi=hi;hi=t;} if(lo>mi){int t=lo;lo=mi;mi=t;}
            if(cyc(lo,mi,hi)!=target[lo][mi][hi]) ok=0; }
        if(ok) dfs(p+1); arc[i][j]=-1; } }
static int cmp(const void*a,const void*b){unsigned long long x=*(unsigned long long*)a,y=*(unsigned long long*)b; return x<y?-1:x>y;}
int main(int argc,char**argv){ N=atoi(argv[1]); char line[256]; long idx=0; np=0; for(int i=0;i<N;i++)for(int j=i+1;j<N;j++){pairs[np][0]=i;pairs[np][1]=j;np++;}
    long bad=0, tot=0;
    while(fgets(line,sizeof line,stdin)){ int L=strlen(line); while(L&&(line[L-1]=='\n'||line[L-1]=='\r')) line[--L]=0; if(L!=np) continue; idx++;
        int p=0; for(int i=0;i<N;i++)for(int j=i+1;j<N;j++) arc[i][j]=(line[p++]=='1');
        for(int a=0;a<N;a++)for(int b=a+1;b<N;b++)for(int c=b+1;c<N;c++) target[a][b][c]=cyc(a,b,c);
        for(int i=0;i<N;i++)for(int j=i+1;j<N;j++) arc[i][j]=-1;
        cnt=0; ns=0; dfs(0); qsort(scores,ns,sizeof(unsigned long long),cmp); int d=0; for(int i=0;i<ns;i++) if(i==0||scores[i]!=scores[i-1]) d++;
        tot++; if(d!=cnt){ bad++; printf("class %ld: C3-class size %lld, distinct score functions %d\n", idx, cnt, d); } }
    printf("n=%d classes=%ld classes where (C3,score function) does NOT determine T: %ld\n", N, tot, bad); return 0; }
