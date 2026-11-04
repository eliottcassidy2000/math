#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;
// Exact physical-subset DP compressed only by independent rotations of each C3.
// No occupied-run normalization and no truncation in child run multiplicity.
int physical(const vector<int>& q,int a,int c) {
 const int k=q.size(), n=k+2, N=1<<(2*k+2); vector<int> incoming(n); vector<int> step(n);
 for(int i=0;i<k;i++) { step[i]=1<<(2*i); for(int j=0;j<k;j++) if(q[j]>>i&1) incoming[i]|=1<<j; if(!(a>>i&1)) incoming[i]|=1<<k; if(!(c>>i&1)) incoming[i]|=1<<(k+1); }
 step[k]=1<<(2*k); step[k+1]=1<<(2*k+1);incoming[k]=a;incoming[k+1]=c|(1<<k);
 const int stride=2*k+2; vector<uint8_t> d((size_t)N*stride);
 for(int s=1;s<N;s++) {
   auto *out=&d[(size_t)s*stride];
   for(int i=0;i<n;i++) {
     int occ=i<k?(s>>(2*i))&3:(s>>(2*k+i-k))&1;if(!occ)continue;
     int prev=s-step[i];const auto* in=&d[(size_t)prev*stride];int ext=prev==0;
     int mask=incoming[i];while(mask){int j=__builtin_ctz((unsigned)mask);mask&=mask-1;ext+=in[j];}
     if(i>=k||occ==1)out[i]=ext&7;
     else if(occ==2){out[i]=(2*ext+in[i])&7;out[n+i]=(ext+in[i])&7;}
     else out[i]=(3*(ext+in[n+i]))&7;
   }
 }
 int answer=0;for(int i=0;i<n;i++)answer+=d[(size_t)(N-1)*stride+i];return answer&7;
}
int main(){int cases;cin>>cases;for(int z=0;z<cases;z++){int id,k,a,c;cin>>id>>k>>a>>c;vector<int>q(k),op(k);for(int &r:q)cin>>r;for(int i=0;i<k;i++)op[i]=((1<<k)-1)^q[i]^(1<<i);auto t=chrono::steady_clock::now();int s=physical(q,a,c),v=physical(op,a,c);cout<<id<<' '<<s<<' '<<v<<'\n';cout.flush();cerr<<z+1<<' '<<id<<' '<<chrono::duration<double>(chrono::steady_clock::now()-t).count()<<" s\n";}}
