#include <cassert>
#include <iostream>
#include <vector>
#include <tuple>
using namespace std;
int n;vector<pair<int,int>> E;
int run(int m,bool *sc=nullptr,int a0=-1,int b0=-1,int a1=-1,int b1=-1){
 int R[8]={};for(int i=0;i<n-1;i++)R[i]|=1<<(i+1);
 for(int e=0;e<(int)E.size();e++){auto[a,b]=E[e];if(m>>e&1)R[b]|=1<<a;else R[a]|=1<<b;}
 if(a0>=0){R[a0]^=1<<b0;R[b0]^=1<<a0;}
 if(a1>=0){R[a1]^=1<<b1;R[b1]^=1<<a1;}
 if(sc){int S=1;for(int k=0;k<n;k++)for(int x=0;x<n;x++)if(S>>x&1)S|=R[x];
 *sc=S==(1<<n)-1;S=1;for(int k=0;k<n;k++)for(int x=0;x<n;x++)if(S&R[x])S|=1<<x;
 *sc&=S==(1<<n)-1;}
 int dp[256][8]={}; for(int i=0;i<n;i++)dp[1<<i][i]=1;
 for(int M=1;M<(1<<n);M++)for(int i=0;i<n;i++)if(dp[M][i])
   for(int a=R[i]&~M;a;a&=a-1){int j=__builtin_ctz(a);dp[M|1<<j][j]+=dp[M][i];}
 int H=0;for(int i=0;i<n;i++)H+=dp[(1<<n)-1][i];return H;
}



void set_order(int a){n=a;E.clear();for(int i=0;i<n;i++)for(int j=i+2;j<n;j++)E.emplace_back(i,j);}
int main(){
 set_order(7);
 vector<vector<int>> C;
 int to59=0,through59to53=0;
 for(int m=0;m<1<<(int)E.size();m++)if(run(m)==55){
    for(int a=0;a<n;a++)for(int b=a+1;b<n;b++){
       int h=run(m,nullptr,a,b);
       assert(h!=51 && h!=53 && h!=57);
       if(h!=59)continue;
       to59++;
       int exits=0;
       for(int c=0;c<n;c++)for(int d=c+1;d<n;d++){
          int t=run(m,nullptr,a,b,c,d);
          assert(t!=57);through59to53+=(t==53); exits+=(t==53);
       }
       assert(exits==1);
    }
    vector<int> R(n);for(int i=0;i<n-1;i++)R[i]|=1<<(i+1);
    for(int e=0;e<(int)E.size();e++){auto[a,b]=E[e];if(m>>e&1)R[b]|=1<<a;else R[a]|=1<<b;}
    C.push_back(R);
 }
 assert(C.size()==220);
 assert(run(81)==55&&run(85)==49&&run(597)==53&&run(85,nullptr,2,3)==57);
 assert(run(113)==59&&run(117)==53);
 // Same marked spine, even for the closer common detour.
 assert(run(628)==55&&run(629)==49&&run(117)==53&&run(565)==57);
 cout<<"Core first arcs to 59: "<<to59<<"; two-edit routes through 59 to 53: "<<through59to53<<endl;
 set_order(8);
 for(auto &R:C) for(int side=0;side<2;side++){
    int m=0;
    for(int i=0;i<(int)E.size();i++){
       auto[a,b]=E[i];bool back=side?(b!=7&&(R[b]>>a&1)):(a>0&&(R[b-1]>>(a-1)&1));
       if(back)m|=1<<i;
    }
    assert(run(m)==55);
    for(int u=0;u<7;u++){
       int a=side?u:0,b=side?7:u+1;
       assert(run(m,nullptr,a,b)!=59);
    }
 }
 cout<<"No padded cross arc to 59. Fixed-spine common setup: 628 → 629 → 117 or 565."<<endl;
}
