#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>
using I=__int128;using U=unsigned __int128;using ull=uint64_t;
using namespace std;
I P(int k){return I(1)<<k;}
I f(I x){I y=0;while(x){I t=x>>1;if(!(x&1))y+=t+(t&1);x=t;}return y;}
string show(U i){if(!i)return "0";string s;while(i){s+=char('0'+i%10);i/=10;}reverse(s.begin(),s.end());return s;}
long long inv(long long a,long long m){
 long long u=0,v=1,z=m;
 while(a){long long q=z/a,t=z%a;z=a;a=t;t=u-q*v;u=v;v=t;}
 assert(z==1);return (u%m+m)%m;
}
I mod(I a,I b){a%=b;return a<0?a+b:a;}
int main(int argc,char**argv){
 int width=argc>1?stoi(argv[1]):54;
 int errors=argc>2?stoi(argv[2]):2;
 assert(width>0&&width<61&&errors>=0&&errors<=3);
 const I target=P(121)-P(60)-1;
 const I B=P(width),full=B-1;
 vector<ull> sparse{0};
 if(errors>=1)for(int i=0;i<width;++i)sparse.push_back(1ULL<<i);
 if(errors>=2)for(int i=0;i<width;++i)for(int j=0;j<i;++j)sparse.push_back(1ULL<<i|1ULL<<j);
 if(errors>=3)for(int i=0;i<width;++i)for(int j=0;j<i;++j)for(int k=0;k<j;++k)sparse.push_back(1ULL<<i|1ULL<<j|1ULL<<k);
 struct Rec{I r,rv,fr;};
 vector<Rec> v;
 for(ull e:sparse)for(I r:{I(e),full^e}){
  I rev=0;for(int i=0;i<width;++i)rev=2*rev+(r>>i&1);
  v.push_back({r,rev,f(r)});
 }
 int n;long long c,h,k0,k1;string mask;ull tests=0, pairs=0;int bank=0;
 while(cin>>n>>c>>h>>k0>>k1>>mask){
  ++bank;
  for(auto [r,rev,fr]:v){
   ++pairs;
   I slope=c*B+h*(full-rev);
   long long a=slope%h,g=gcd(a,h),m=h/g,j=m>1?inv(a/g,m):0;
   for(int i=0;i<2;++i){
    I residual=target-(r&1?k1:k0)-c*r-h*(fr+(i&&r<B/2));
    if(residual<0||residual%g)continue;
    I left=(residual+slope+h-1)/(slope+h),right=residual/slope;
    // Width means an actual low-width suffix; leading-zero padded addresses are excluded.
    left=max(left,I(1));
    if(left>right)continue;
    I res=m==1?0:(residual/g%m)*j%m;
    if((res&1)!=i)res+=m; // m odd for receiver charts
    assert(m&1);
    long long step=2*m;
    for(I z=left+mod(res-left,step);z<=right;z+=step){
     ++tests;
     if(slope*z+h*f(z)!=residual)continue;
     I y=B*z+r;
     assert((y&1?k1:k0)+c*y+h*f(y)==target);
     cout<<"HIT "<<mask<<" "<<show(y)<<"\n";
    }
   }
  }
 }
 cerr<<"width="<<width<<", errors="<<errors<<", bank="<<bank<<", suffix/receiver pairs="<<pairs<<", exact prefix tests="<<tests<<"\n";
}
