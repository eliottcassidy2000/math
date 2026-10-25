#include <cstdint>
#include <iostream>
#include <vector>
#include <stdexcept>
// Fixed start 0, remaining-vertex mask. Bound n<=19 is intentional.
std::vector<uint64_t> count(int n, bool backward) {
 if(n<3||n>19||n%4!=3) throw std::runtime_error("outside certified range");
 int m=n-1; uint32_t size=1u<<m;
 std::vector<bool> qr(n);for(int x=1;x<n;x++)qr[x*x%n]=true;
 std::vector<uint32_t> out(m),in(m);
 for(int v=1;v<n;v++)for(int w=1;w<n;w++)if(qr[(w-v+n)%n]) {out[v-1]|=1u<<(w-1);in[w-1]|=1u<<(v-1);}
 std::vector<uint64_t> dp(uint64_t(size)*m);
 for(int v=1;v<n;v++)if(qr[v])dp[uint64_t(1u<<(v-1))*m+v-1]=1;
 for(uint32_t mask=1;mask<size;mask++) {
   uint32_t vs=mask;while(vs){int v=__builtin_ctz(vs);vs&=vs-1;
    if(backward) {
      uint32_t prev=mask^(1u<<v);if(!prev)continue;
      uint32_t us=prev&in[v]; uint64_t total=0;
      while(us){int u=__builtin_ctz(us);us&=us-1;total+=dp[uint64_t(prev)*m+u];}
      dp[uint64_t(mask)*m+v]=total;
    } else {
      uint64_t val=dp[uint64_t(mask)*m+v]; if(!val)continue;
      uint32_t ws=out[v]&~mask;
      while(ws){int w=__builtin_ctz(ws);ws&=ws-1;dp[uint64_t(mask|(1u<<w))*m+w]+=val;}
    }
   }
 }
 std::vector<uint64_t> ans(n);for(int v=1;v<n;v++)ans[v]=dp[uint64_t(size-1)*m+v-1];return ans;
}
int main(){for(int n:{3,7,11,19}){auto a=count(n,false),b=count(n,true);if(a!=b)throw std::runtime_error("recurrence mismatch");std::cout<<n;for(auto x:a)std::cout<<" "<<x;std::cout<<"\n";}}
