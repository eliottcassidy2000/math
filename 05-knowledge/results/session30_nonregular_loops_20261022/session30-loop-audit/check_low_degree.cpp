#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>
using namespace std;
using I=long long;
int main(){
 for(int n=1;n<=7;n++){
  int ne=n*(n-1)/2; uint64_t universe=1ULL<<ne;
  uint64_t cases=n<=6?universe:128, ar=0; I min4=0,max4=0;
  vector<int> p(n); iota(p.begin(),p.end(),0); vector<vector<int>> perms;
  do{perms.push_back(p);}while(next_permutation(p.begin(),p.end()));
  I fact[9]={1};for(int j=1;j<=8;j++)fact[j]=fact[j-1]*j;
  uint64_t rng=0x23dd734f;
  for(uint64_t c=0;c<cases;c++){
   uint64_t mask=c;if(n==7){rng^=rng<<13;rng^=rng>>7;rng^=rng<<17;mask=rng%universe;}
   I s[8][8]={},d[8]={},sd[8]={};int bit=0;
   for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){s[i][j]=(mask>>bit++&1)?1:-1;s[j][i]=-s[i][j];}
   for(int i=0;i<n;i++)for(int j=0;j<n;j++)d[i]+=s[i][j];
   for(int i=0;i<n;i++)for(int j=0;j<n;j++)sd[i]+=s[i][j]*d[j];
   I ds=0,d4=0,ds2d=0,d2sd=0;
   bool almost=n%2==0;
   for(int i=0;i<n;i++){ds+=d[i]*d[i];d4+=d[i]*d[i]*d[i]*d[i];ds2d-=sd[i]*sd[i];d2sd+=d[i]*d[i]*sd[i];almost &= abs(d[i])==1;}
   I sum[5]={};
   for(const auto&q:perms){I poly[5]={1,0,0,0,0};for(int j=0;j<n-1;j++)for(int l=min(4,j+1);l>=1;l--)poly[l]+=s[q[j]][q[j+1]]*poly[l-1];for(int l=0;l<5;l++)sum[l]+=poly[l];}
   I N[5]={};for(int l=0;l<5;l++)if(l<n){if(sum[l]%fact[n-l])return 4;N[l]=sum[l]/fact[n-l];}
   I predicted2=I(n)*(n-1)-ds;
   I predicted4=I(n)*(n-1)*(15-7*n)+2*(8*n-13)*ds+2*ds2d+4*d2sd-d4;
   I observed4=2*N[4]-N[2]*N[2];
   if(N[2]!=predicted2 || observed4!=predicted4 || N[1] || N[3]){cerr<<"FAIL n="<<n<<" mask="<<mask<<" N2="<<N[2]<<" pred2="<<predicted2<<" z4numer="<<observed4<<" pred4="<<predicted4<<"\n";return 1;}
   if(almost){ar++;I simplified=I(n)*(n-1)*(15-7*n)+2*ds2d+I(n)*(16*n-31);if(simplified!=observed4)return 2;}
   if(c==0)min4=max4=observed4;else{min4=min(min4,observed4);max4=max(max4,observed4);}
  }
  cout<<"n="<<n<<" cases="<<cases<<" scope="<<(n<=6?"ALL_LABELLED":"SEEDED_SAMPLE")<<" almost_regular="<<ar<<" numerator4_range="<<min4<<","<<max4<<" PASS\n";
 }
}
