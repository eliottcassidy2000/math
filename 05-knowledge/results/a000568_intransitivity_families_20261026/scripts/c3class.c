/* c3class.c -- for every tournament class (gentourng input), count the tournaments on the SAME labelled
 * vertex set with the SAME set of cyclic triples (C3-structure), by DFS over arcs with triple checks.
 * Prints per class: idx, count. Also prints the tree prediction prod(k!)*2^prime computed independently
 * (from the table produced by tourn_enum) is compared downstream in Python.
 * Build: cc -O2 c3class.c -o c3class ; use: gentourng n | ./c3class n */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int N; static int target[12][12][12]; /* 1 if triple cyclic */
static int arc[12][12]; /* arc[i][j]=1 if i->j, for i<j assigned; -1 unknown */
static int pairs[66][2], np;
static long long cnt;
static int cyc(int a,int b,int c){ /* uses arc[] for all three pairs (assigned) */
    int ab = a<b? arc[a][b] : !arc[b][a]; int bc = b<c? arc[b][c] : !arc[c][b]; int ca = c<a? arc[c][a] : !arc[a][c];
    return (ab&&bc&&ca)||(!ab&&!bc&&!ca);
}
static void dfs(int p){
    if(p==np){cnt++;return;}
    int i=pairs[p][0], j=pairs[p][1];
    for(int d=0; d<2; d++){
        arc[i][j]=d;
        int ok=1;
        /* check all triples (i,j,k) with pair (i,k),(j,k) already assigned: pairs are in lexicographic order (i<j),
           so (i,k) and (j,k) assigned iff index < p. We check all k != i,j whose both pairs assigned. */
        for(int k=0;k<N&&ok;k++){ if(k==i||k==j) continue;
            int a=i<k?i:k, b=i<k?k:i; int c=j<k?j:k, e=j<k?k:j;
            if(arc[a][b]<0||arc[c][e]<0) continue;
            int x=i,y=j,z=k; int lo=x,mi=y,hi=z; /* sort */
            if(lo>mi){int t=lo;lo=mi;mi=t;} if(mi>hi){int t=mi;mi=hi;hi=t;} if(lo>mi){int t=lo;lo=mi;mi=t;}
            if(cyc(lo,mi,hi)!=target[lo][mi][hi]) ok=0;
        }
        if(ok) dfs(p+1);
        arc[i][j]=-1;
    }
}
int main(int argc,char**argv){ N=atoi(argv[1]); char line[256]; long idx=0;
    np=0; for(int i=0;i<N;i++)for(int j=i+1;j<N;j++){pairs[np][0]=i;pairs[np][1]=j;np++;}
    while(fgets(line,sizeof line,stdin)){ int L=strlen(line); while(L&&(line[L-1]=='\n'||line[L-1]=='\r')) line[--L]=0; if(L!=np) continue; idx++;
        int p=0; for(int i=0;i<N;i++)for(int j=i+1;j<N;j++) arc[i][j]=(line[p++]=='1');
        for(int a=0;a<N;a++)for(int b=a+1;b<N;b++)for(int c=b+1;c<N;c++) target[a][b][c]=cyc(a,b,c);
        for(int i=0;i<N;i++)for(int j=i+1;j<N;j++) arc[i][j]=-1;
        cnt=0; dfs(0); printf("%ld %lld\n", idx, cnt);
    }
    return 0; }
