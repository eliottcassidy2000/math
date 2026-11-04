// Independent backward completion counters. No producer code is linked.
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;
using U=uint64_t;

// F(state, first) counts suffix words with exactly the remaining multiplicities.
// Removing the first symbol strictly decreases the mixed-radix state.
vector<U> starts(const vector<unsigned>& arc,const vector<int>& cap){
 int n=arc.size(), states=1; vector<int> weight(n);
 for(int i=0;i<n;++i){weight[i]=states;states*=cap[i]+1;}
 vector<unsigned> present(states,0);
 for(int s=1;s<states;++s) for(int j=0;j<n;++j)
  if((s/weight[j])%(cap[j]+1)) present[s]|=1u<<j;
 vector<U> f(size_t(states)*n,0);
 for(int s=1;s<states;++s){
  unsigned possible=present[s];
  while(possible){int i=__builtin_ctz(possible);possible&=possible-1;
   int rest=s-weight[i]; U ans=0;
   if(rest==0) ans=1;
   else {unsigned next=arc[i]&present[rest];while(next){int j=__builtin_ctz(next);next&=next-1;ans+=f[size_t(rest)*n+j];}}
   f[size_t(s)*n+i]=ans;
  }
 }
 return vector<U>(f.end()-n,f.end());
}

// Literal distinct-vertex suffix recurrence; optional root closure counts cycles.
U physical(const vector<unsigned>& arc,int root){
 int n=arc.size(); unsigned end=1u<<n; vector<U> f(size_t(end)*n,0);
 for(unsigned m=1;m<end;++m){
  unsigned initial=m;
  while(initial){int i=__builtin_ctz(initial);initial&=initial-1;
   unsigned rest=m^(1u<<i);U ans=0;
   if(!rest) ans=(root<0 || ((arc[i]>>root)&1));
   else {unsigned next=arc[i]&rest;while(next){int j=__builtin_ctz(next);next&=next-1;ans+=f[size_t(rest)*n+j];}}
   f[size_t(m)*n+i]=ans;
  }
 }
 U ans=0;
 if(root<0){for(int i=0;i<n;++i)ans+=f[size_t(end-1)*n+i];}
 else {unsigned rem=(end-1)^(1u<<root),next=arc[root]&rem;while(next){int j=__builtin_ctz(next);next&=next-1;ans+=f[size_t(rem)*n+j];}}
 return ans;
}
int main(){ios::sync_with_stdio(false);cin.tie(nullptr);int jobs;cin>>jobs;
 while(jobs--){char mode;int n;cin>>mode>>n;vector<unsigned> arc(n);for(auto&x:arc)cin>>x;
  if(mode=='R'||mode=='A'){vector<int> cap(n);for(auto&x:cap)cin>>x;auto values=starts(arc,cap);
   if(mode=='A'){for(int i=0;i<n;++i){if(i)cout<<' ';cout<<values[i];}cout<<'\n';}
   else {U sum=0;for(auto x:values)sum+=x;cout<<sum<<'\n';}
  }else if(mode=='H'){cout<<physical(arc,-1)<<'\n';}
  else if(mode=='Z'){int root;cin>>root;cout<<physical(arc,root)<<'\n';}
  else return 2;
 }
}
