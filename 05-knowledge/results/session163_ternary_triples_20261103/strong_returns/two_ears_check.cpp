#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cstdlib>
using namespace std;
void need(bool p,const char*m){if(!p){fprintf(stderr,"FAIL: %s\n",m);exit(1);}}
static long dp[16384][14];
long H(const int*a,int n){int N=1<<n;memset(dp,0,N*sizeof(dp[0]));for(int i=0;i<n;i++)dp[1<<i][i]=1;for(int m=1;m<N;m++)for(int t=m;t;t&=t-1){int i=__builtin_ctz((unsigned)t);long x=dp[m][i];if(x)for(int z=a[i]&~m;z;z&=z-1){int j=__builtin_ctz((unsigned)z);dp[m|1<<j][j]+=x;}}long s=0;for(int i=0;i<n;i++)s+=dp[N-1][i];return s;}
void add(int*a,int n,int port){a[n]=1<<port;for(int i=0;i<n;i++)if(i!=port)a[i]|=1<<n;}
bool strong(int*a,int n){for(int i=0;i<n;i++){int z=1<<i,old=0;while(old!=z){old=z;for(int j=0;j<n;j++)if(z>>j&1)z|=a[j];}if(z!=(1<<n)-1)return false;}return true;}
void attach(int*a,const int*x,int n,int s){for(int i=0;i<n;i++)a[i]=x[i]|((s>>i&1)?0:1<<n);a[n]=s;}

int main(){const int src[2][10]={{1008,1017,1011,997,776,912,816,592,128,256},{996,993,994,999,975,912,800,576,128,256}};
for(int id=0;id<2;id++){
 int x[14],z[14],good[1024],gc=0,cap=id?877:859;long targets[3]={id?873:855,id?875:857,id?877:859};
 for(int s=0;s<1024;s++){attach(x,src[id],10,s);if(H(x,11)<=cap)good[gc++]=s;}
 need(gc==(id?46:32),"surviving one-ear masks");
 long cnt=0,hits=0,shits=0,htarget=0; long scalar[3]={},st[3]={},observed[3]={};
 for(int i=0;i<gc;i++)for(int j=0;j<gc;j++){
  attach(x,src[id],10,good[i]);x[10]|=1<<11;for(int u=0;u<10;u++)x[u]|=((good[j]>>u&1)?0:1<<11);x[11]=good[j];
  long h=H(x,12);cnt++;if(h!=targets[0]&&h!=targets[1]&&h!=targets[2])continue;htarget++;int tid=h==targets[0]?0:(h==targets[1]?1:2);scalar[tid]++;if(!strong(x,12))continue;shits++;st[tid]++;
  for(int p=0;p<12;p++){copy(x,x+12,z);add(z,12,p);if(H(z,13)!=1593)continue;add(z,13,12);if(H(z,14)!=2943)continue;
    printf("MATCH source=%d target=%ld port=%d A=%d B=%d rows",id,h,p,good[i],good[j]);for(int u=0;u<12;u++)printf(" %d",x[u]);puts("");hits++;observed[tid]++;
  }
 }
 need(cnt==(id?2116:1024),"pair count");
 need(scalar[0]==(id?6:7) && scalar[1]==0 && scalar[2]==0,"scalar target census");
 need(st[0]==scalar[0] && st[1]==0 && st[2]==0,"strong target census");
 need(observed[0]==(id?3:5) && observed[1]==0 && observed[2]==0,"marked observer census");
 printf("scalar [%ld,%ld,%ld], strong [%ld,%ld,%ld], full states [%ld,%ld,%ld]\n",scalar[0],scalar[1],scalar[2],st[0],st[1],st[2],observed[0],observed[1],observed[2]);
 printf("id %d, retained singles %d, pairs %ld, target scalar %ld, strong scalar %ld, full observer %ld\n",id,gc,cnt,htarget,shits,hits);fflush(stdout);
}
 puts("PASS");
}
