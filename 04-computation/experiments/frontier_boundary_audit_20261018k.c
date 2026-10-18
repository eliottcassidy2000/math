/* Independent every-integer audit of the fixed frontier source. No route bank,
   no event-order theorem, and no reuse of the Python shell construction. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
int main(int argc,char **argv){
 assert(argc==3);int power=atoi(argv[2]);assert(power>=3 && power<=30);
 int64_t lim=1LL<<power;int v[4096],nv=0,maxv=0;char line[256];
 FILE *f=fopen(argv[1],"r");assert(f);
 while(fgets(line,sizeof(line),f)){int x;if(line[0]!='#' && sscanf(line,"%d",&x)==1){assert(nv<4096);v[nv++]=x;if(x>maxv)maxv=x;}}
 fclose(f);assert(nv==627);
 int *index=malloc((maxv+1)*sizeof(int));assert(index);
 for(int i=0;i<=maxv;i++)index[i]=-1;
 for(int i=0;i<nv;i++)index[v[i]]=i;
 int c[4096][3]={{0}};int64_t D=0,mx=INT64_MIN,mn=INT64_MAX,mxn=0,mnn=0,updates=0;
 printf("{\"status\":\"FINITE-EXACT\",\"limit_exclusive\":%lld,\"dyadics\":[",(long long)lim);int sep=0;
 for(int64_t N=1;N<lim;N++){
   for(int event=0;event<3;event++){
     int64_t m=event==0?4*N-3:event==1?4*N-1:2*N-1,x=m;
     while(x%8==5)x=(x-1)/4;
     if(x>maxv || index[x]<0)continue;
     int i=index[x];
     for(int r=0;r<3;r++){
       int residue=(m+r)%3,d=event==2?(residue==1?-2:1):(residue==2?2:-1);
       D+=(r==0?2:-1)*((int64_t)2*c[i][r]*d+d*d);c[i][r]+=d;
       assert(c[i][r]>=-3 && c[i][r]<=3);
     }
     updates++;
   }
   if(D>mx){mx=D;mxn=N;}if(D<mn){mn=D;mnn=N;}
   if((N&(N-1))==0)printf("%s[%lld,%lld]",sep++?",":"",(long long)N,(long long)D);
 }
 printf("],\"maximum\":%lld,\"maximum_N\":%lld,\"minimum\":%lld,\"minimum_N\":%lld,\"updates\":%lld}\n",(long long)mx,(long long)mxn,(long long)mn,(long long)mnn,(long long)updates);
 free(index);return 0;
}
