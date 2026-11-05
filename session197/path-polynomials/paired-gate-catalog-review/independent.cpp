#include <algorithm>
#include <array>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>
using namespace std; using ll=long long; using Graph=vector<int>;
void need(bool t,const char *s){if(!t)throw runtime_error(s);}
int pc(int x){return __builtin_popcount((unsigned)x);}
Graph canonical(const Graph&a){int n=a.size();vector<vector<int>>cells(n);for(int v=0;v<n;v++)cells[pc(a[v])].push_back(v);Graph p,best;
 function<void(int)> rec=[&](int d){if(d==n){Graph rows(n);for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(a[p[i]]>>p[j]&1)rows[i]|=1<<j;if(best.empty()||rows<best)best=rows;return;}auto cell=cells[d];do{p.insert(p.end(),cell.begin(),cell.end());rec(d+1);p.resize(p.size()-cell.size());}while(next_permutation(cell.begin(),cell.end()));};rec(0);return best;}
ll aut(const Graph&a){int n=a.size();vector<vector<int>>cells(n);for(int v=0;v<n;v++)cells[pc(a[v])].push_back(v);Graph p(n,-1);ll count=0;
 function<void(int)> rec=[&](int d){if(d==n){for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(((a[i]>>j)&1)!=((a[p[i]]>>p[j])&1))return;count++;return;}auto c=cells[d],q=c;do{for(int i=0;i<(int)c.size();i++)p[c[i]]=q[i];rec(d+1);}while(next_permutation(q.begin(),q.end()));};rec(0);return count;}
array<ll,3> covers(const Graph&a){int n=a.size();vector<ll>h(1<<n);function<void(int,int)> path=[&](int S,int v){h[S]++;for(int bits=a[v]&~S;bits;bits&=bits-1){int b=bits&-bits;path(S|b,__builtin_ctz((unsigned)b));}};for(int i=0;i<n;i++)path(1<<i,i);
 array<ll,3> out{};vector<int>blocks;function<void(int)> partition=[&](int v){if(v==n){ll c=1;for(int S:blocks)c*=h[S];out[blocks.size()-1]+=c;return;}for(int j=0;j<(int)blocks.size();j++){blocks[j]|=1<<v;partition(v+1);blocks[j]^=1<<v;}if(blocks.size()<3){blocks.push_back(1<<v);partition(v+1);blocks.pop_back();}};partition(0);return out;}
ll scalar(array<ll,3>x){ll h=x[0],u=x[1],v=x[2];return 5*h*h+12*h*u+28*u*u+48*u*v+72*v*v;}
ll powi(ll a,int k){ll b=1;while(k--)b*=a;return b;}
ll factor(ll x){if(x<2)return 0;for(ll q=2;q<=x/q;q++)if(x%q==0)return q;return 1;}
Graph pair_graph(const Graph&a){int n=a.size();Graph g(2*n+2);for(int i=0;i<n;i++){g[i]=a[i]|(((1<<n)-1)<<n);g[n+i]=(((1<<n)-1)^a[i]^(1<<i))<<n;g[n+i]|=3<<(2*n);}g[2*n]=((1<<n)-1)|(1<<(2*n+1));g[2*n+1]=(1<<n)-1;return g;}
ll walk_ie(const Graph&a){int n=a.size(),all=(1<<n)-1;__int128 total=0;for(int S=1;S<=all;S++){__int128 d[16]={},q[16]={};for(int v=0;v<n;v++)if(S>>v&1)d[v]=1;for(int t=1;t<n;t++){fill(q,q+n,0);for(int v=0;v<n;v++)if(d[v])for(int bits=a[v]&S;bits;bits&=bits-1)q[__builtin_ctz((unsigned)bits)]+=d[v];copy(q,q+n,d);}__int128 count=accumulate(d,d+n,(__int128)0);total+=((n-pc(S))%2?-count:count);}need(total>=0&&total<((__int128)1<<63),"walk result range");return (ll)total;}
int main(int argc,char**argv){need(argc==2,"input CSV");set<Graph> bank{{0}};set<ll>image;vector<int>classes;vector<ll>masses;ll fac=1;int literal=0;for(int n=1;n<=7;n++){if(n>1){set<Graph>next;for(auto&a:bank)for(int mask=0;mask<(1<<(n-1));mask++){Graph b=a;b.push_back(mask);for(int i=0;i<n-1;i++)if(!(mask>>i&1))b[i]|=1<<(n-1);next.insert(canonical(b));}bank=move(next);}fac*=n;ll mass=0;for(auto&a:bank){ll au=aut(a);need(fac%au==0,"orbit divisor");mass+=fac/au;ll H=scalar(covers(a));if(H<73529109)image.insert(H);literal++;}need(mass==(1LL<<(n*(n-1)/2)),"complete orbit mass");classes.push_back(bank.size());masses.push_back(mass);}
 need(classes==vector<int>({1,1,2,4,12,56,456}),"unlabelled inventory");ll prev=0;for(int n=1;n<=12;n++){Graph a(n);for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)a[i]|=1<<j;ll value=scalar(covers(a));ll formula=2*powi(9,n)-8*powi(6,n)+13*powi(4,n)+4*powi(3,n)-22*powi(2,n)+15;need(value==formula&&value>prev,"physical Phi controls");prev=value;if(n==8)need(value==73529109,"exclusive cutoff");}
 ifstream in(argv[1]);string line;set<ll>given;int witnesses=0,outputs=0;vector<ll>primes,sevens;while(getline(in,line)){replace(line.begin(),line.end(),',',' ');istringstream ss(line);ll H,x,y,z,f,f7;int n;ss>>H>>x>>y>>z>>f>>f7>>n;Graph a(n);for(int&i:a)ss>>i;need(bool(ss),"CSV parse");for(int i=0;i<n;i++){need(!(a[i]>>i&1),"loop");for(int j=i+1;j<n;j++)need(((a[i]>>j)&1)+((a[j]>>i)&1)==1,"tournament");}auto c=covers(a);need(c==array<ll,3>{x,y,z}&&scalar(c)==H,"literal witness covers");need(factor(H)==f&&(H%7?0:factor(H/7))==f7,"prime fields");given.insert(H);witnesses++;if(f==1)primes.push_back(H);if(f7==1)sevens.push_back(H);if(f==1||f7==1){Graph g=pair_graph(a);vector<int>tau;for(int i=0;i<n;i++)tau.push_back(n+i);for(int i=0;i<n;i++)tau.push_back(i);tau.push_back(2*n+1);tau.push_back(2*n);for(int i=0;i<2*n+2;i++){int S=1<<i,last=0;while(S!=last){last=S;for(int j=0;j<2*n+2;j++)if(S>>j&1)S|=g[j];}need(S==(1<<(2*n+2))-1,"strong paired graph");for(int j=0;j<2*n+2;j++)need(((g[i]>>j)&1)==((g[tau[j]]>>tau[i])&1),"anti-involution");}need(walk_ie(g)==H,"independent walk inclusion exclusion");outputs++;}}
 need(given==image&&image.size()==119,"complete independent image equality");need(primes.size()==8&&sevens==vector<ll>({37050461,46937093,56228501}),"prime image");
 cout<<"{\n  \"status\":\"PASS\",\n  \"classes\":[";for(int i=0;i<7;i++)cout<<(i?",":"")<<classes[i];cout<<"],\n  \"labelled_orbit_masses\":[";for(int i=0;i<7;i++)cout<<(i?",":"")<<masses[i];cout<<"],\n  \"literal_class_cover_counts\":"<<literal<<",\n  \"input_witnesses\":"<<witnesses<<",\n  \"direct_walk_ie_outputs\":"<<outputs<<",\n  \"values\":[";bool first=true;for(ll H:image){cout<<(first?"":",")<<H;first=false;}cout<<"],\n  \"primes\":[";for(int i=0;i<(int)primes.size();i++)cout<<(i?",":"")<<primes[i];cout<<"],\n  \"seven_primes\":[";for(int i=0;i<(int)sevens.size();i++)cout<<(i?",":"")<<sevens[i];cout<<"]\n}\n";
}
