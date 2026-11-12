/* one-vertex extensions of the order-8 classes: print the order-9 tournaments whose H is in the target list
   (argv[1..]) and which are strong. stdin: 8 out-masks per line. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int out[9]; static long long dp[512][9];
static long long H9(void){ memset(dp,0,sizeof(dp)); for(int v=0;v<9;v++) dp[1<<v][v]=1;
    for(int S=1;S<512;S++) for(int v=0;v<9;v++){ long long c=dp[S][v]; if(!c) continue; int t=out[v]&~S; while(t){int b=t&-t; dp[S|b][__builtin_ctz(b)]+=c; t^=b;} }
    long long h=0; for(int v=0;v<9;v++) h+=dp[511][v]; return h; }
static int strong(void){ int r[9]; for(int i=0;i<9;i++) r[i]=out[i]|(1<<i); for(int k=0;k<9;k++) for(int i=0;i<9;i++) if((r[i]>>k)&1) r[i]|=r[k]; for(int i=0;i<9;i++) if(r[i]!=511) return 0; return 1; }
int main(int argc,char**argv){ long long tg[64]; int nt=argc-1; for(int i=0;i<nt;i++) tg[i]=atoll(argv[i+1]);
    char line[256]; int base[8]; long cnt=0;
    while(fgets(line,sizeof(line),stdin)){ char*p=line; int ok=1; for(int i=0;i<8;i++){ char*e; long v=strtol(p,&e,10); if(e==p){ok=0;break;} base[i]=(int)v; p=e; } if(!ok) continue;
        for(int X=0;X<256;X++){ for(int i=0;i<8;i++){ out[i]=base[i]; if(!((X>>i)&1)) out[i]|=1<<8; } out[8]=X;
            long long h=H9(); int hit=0; for(int i=0;i<nt;i++) if(h==tg[i]) hit=1; if(!hit) continue; if(!strong()) continue;
            for(int i=0;i<9;i++) printf("%d ",out[i]); printf("\n"); cnt++; } }
    fprintf(stderr,"order-9 strong tournaments (labelled extensions) with target H: %ld\n",cnt); return 0; }
