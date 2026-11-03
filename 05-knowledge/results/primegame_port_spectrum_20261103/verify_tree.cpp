#include <cstdint>
#include <iostream>
#include <vector>
#include <algorithm>
using I=int64_t;
std::vector<unsigned> a;std::vector<I> memo;int n;
I start(unsigned s,int v){
  I &value=memo[size_t(s)*n+v];if(value>=0)return value;
  unsigned rem=s^(1u<<v);if(!rem)return value=1;
  value=0;unsigned next=a[v]&rem;
  while(next){int w=__builtin_ctz(next);next&=next-1;value+=start(rem,w);}return value;
}
int main(){I target;size_t checked=0;I maxvalue=0;
  while(std::cin>>target>>n){
    if(n<1||n>11)return 2;a.assign(n,0);
    for(int z=1;z<n;z++){int x;std::cin>>x;if(x<0||x>=z)return 3;
      a[z]=1u<<x;for(int i=0;i<z;i++)if(i!=x)a[i]|=1u<<z;
    }
    memo.assign(size_t(1u<<n)*n,-1);I total=0;
    for(int v=0;v<n;v++)total+=start((1u<<n)-1,v);
    if(total!=target||target==7||target==21||target%2!=1){std::cerr<<"failed "<<target<<" "<<total<<"\n";return 4;}
    checked++;maxvalue=std::max(maxvalue,target);
  }
  std::cout<<"{\"verified_witnesses\":"<<checked<<",\"maximum\":"<<maxvalue<<",\"algorithm\":\"top-down starting-endpoint subset recurrence\"}\n";
}
