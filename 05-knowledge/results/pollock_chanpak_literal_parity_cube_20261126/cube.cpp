#include <iostream>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
using namespace std;
int H(int n, uint32_t m, int &cyc, bool flag){
 int out[10]={}; int k=0;
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++)
  if((m>>k)&1)out[j]|=1<<i;else out[i]|=1<<j;
 if(flag){int s=1,t;do{t=s;for(int j=0;j<n;j++)if(s>>j&1)s|=out[j];}while(t!=s);
 if(s!=(1<<n)-1)cyc=0;
 s=1;do{t=s;for(int j=0;j<n;j++)if(out[j]&s)s|=1<<j;}while(t!=s);
 if(s!=(1<<n)-1)cyc=0;}
 int dp[1<<8][8]={};
 for(int i=0;i<n;i++)dp[1<<i][i]=1;
 for(int s=1;s<(1<<n);s++)for(int i=0;i<n;i++)if(dp[s][i]){
  int rem=out[i]&~s;
  while(rem){int j=__builtin_ctz(rem);rem &= rem-1;dp[s|(1<<j)][j]+=dp[s][i];}
 }
 int h=0;for(int i=0;i<n;i++)h+=dp[(1<<n)-1][i]; return h;
}
int main(int argc,char**argv){
 int n=atoi(argv[1]),k=0,a[3];
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++)if(i%2==0&&j==i+1&&i<6)a[i/2]=k;
 int basis[8]={}, mark=0;
 for(int x=1;x<8;x++) for(int j=0;j<3;j++)if(x>>j&1)basis[x]|=1<<a[j];
 mark=basis[7];
 int hist[10001]={};int best=10001, pure=0,par=0,pureStrong=0,pStrong=0,unit=0;
 for(unsigned m=0;m<(1u<<k);m++)if(!(m&mark)){
  int h[8],cyc=1,dp=0;
  for(int x=0;x<8;x++) {h[x]=H(n,m|basis[x],cyc,true);dp+=(__builtin_parity(x)?1:-1)*h[x];}
  if(dp&&abs(dp)<=10000) hist[abs(dp)]++;
  if(abs(dp)==2) unit++;
  bool parity=dp&&h[0]==h[3]&&h[0]==h[5]&&h[0]==h[6]&&h[1]==h[2]&&h[1]==h[4]&&h[1]==h[7];
  if(parity){par++;if(cyc)pStrong++;if(par==1|| (cyc&&pStrong==1)){cout<<"parity "<<m<<" "<<cyc<<" ";for(int z:h)cout<<z<<",";cout<<endl;}}
  if(dp&&abs(dp)<best){best=abs(dp);cout<<"least "<<m<<" "<<cyc<<" ";for(int z:h)cout<<z<<",";cout<<endl;}
  if(dp) for(int s=0;s<8;s++){bool good=h[s]!=h[s^7];for(int x=1;x<7;x++)good&=(h[s]==h[s^x]);
   if(good){pure++;if(cyc)pureStrong++;if(pure==1|| (cyc&&pureStrong==1)){cout<<"pure "<<m<<" "<<s<<" "<<cyc<<" ";for(int z:h)cout<<z<<",";cout<<endl;}}
  }
 }
 cout<<"n "<<n<<" pure "<<pure<<" strong_pure "<<pureStrong<<" parity "<<par<<" strong_parity "<<pStrong<<" unit "<<unit<<" min "<<best<<endl;
}
