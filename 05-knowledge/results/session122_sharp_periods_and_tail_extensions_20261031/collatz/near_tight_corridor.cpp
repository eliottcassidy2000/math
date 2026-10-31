// Independent of producer scheduling: reconstruct all reverse edges with integer slack<=2.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <stdexcept>
#include <limits>
using namespace std;using ll=long long;using Pair=pair<ll,ll>;
struct Node{int z,m,b;vector<pair<int,int>> pred;vector<Pair> pairs;uint64_t paths=0;};
ll power(ll a,ll n,ll q){ll x=1;while(n){if(n&1)x=(__int128)x*a%q;a=(__int128)a*a%q;n>>=1;}return x;}
int main(int argc,char**argv){try{
 if(argc<3)throw runtime_error("E binary [R=2] [budget=2] [endpoint=E]");int E=stoi(argv[1]),R=argc>3?stoi(argv[3]):2,B=argc>4?stoi(argv[4]):2,END=argc>5?stoi(argv[5]):E;
 if(E<2||E>10||R<1||R>E||B<0||B>2||END<E||END>20)throw runtime_error("parameters");
 vector<int>p(E+1,1);for(int i=1;i<=E;i++)p[i]=5*p[i-1];int q=p[E];vector<uint8_t>h(q,E);for(int z=1;z<q;z++){int x=z;h[z]=0;while(x%5==0)x/=5,++h[z];}
 vector<int>base(E+1),width(E+1);ll size=0;for(int e=0;e<=E;e++){base[e]=size;width[e]=4*p[max(R,e)-1];size+=ll(e==E?1:4*q/p[e+1])*width[e];if(size>numeric_limits<int>::max())throw runtime_error("state index overflow");}
 auto index=[&](int z,int m){int e=h[z],j=z?z/p[e]:0;return base[e]+(z?(j-1-j/5)*width[e]:0)+m-1-m/5;};
 ifstream f(argv[2],ios::binary|ios::ate);if(!f||f.tellg()!=8*size)throw runtime_error("binary size");f.seekg(0);vector<int16_t>V(size);array<unsigned char,65536>buf;ll at=0;
 while(at<size){int n=min<ll>(8192,size-at);f.read((char*)buf.data(),8*n);if(!f)throw runtime_error("binary read");for(int i=0;i<n;i++){uint64_t u=0;for(int j=7;j>=0;j--)u=(u<<8)|buf[8*i+j];ll x=(ll)u;if(x<=-30000||x>0)throw runtime_error("potential range");V[at+i]=x;}at+=n;}
 int lo=*min_element(V.begin(),V.end()),hi=*max_element(V.begin(),V.end());
 if(hi+8*E-325>lo-B-1||8*E-325>lo-B-1)throw runtime_error("omitted letter gap not strict enough");
 auto value=[&](int z,int m){return int(V[index(z,m)]);};auto key=[&](int z,int m,int b){return uint64_t(index(z,m))*(B+1)+b;};
 auto ones=[&](ll k){ll s=0;for(int e=1;e<=E;e++)s+=k/(2*p[e-1]);return s;};
 ll inv3=(2LL*q+1)/3; // 5^E is 1 or2 mod3, use generic inverse below.
 inv3=power(3,4LL*q/5-1,q);vector<int>p2(65,1);for(int a=1;a<=64;a++)p2[a]=2LL*p2[a-1]%q;
 unordered_map<uint64_t,vector<Pair>>entries;ll entrytests=0,longtests=0;
 for(int a=2;a<=64;a++){ll inv=power(p2[a],4LL*q/5-1,q);int z=(q+1-2*inv%q)%q,e=h[z],ell=max(R,e),mod=p[ell],T=2*p[ell-1];int rho=3LL*((mod+1)/2)%mod,m=1;
  for(int k=0;k<T;k++){int m2=3LL*m*inv%mod,reward=min<int>(e,h[m2-1]);ll score=8*(ones(k)+reward)-5LL*(k+a),gap=value(z,m2)-score;if(gap<0)throw runtime_error("bad entry inequality");
   if(gap<=B)entries[key(z,m2,gap)].push_back({k+a,k+1});
   if(gap+2<=B)for(int e2=ell;e2<=E;e2++){ll K=k+2*p[e2-1];if(8*(ones(K)+reward)-5*(K+a)!=score-2)throw runtime_error("long entry score");if(3LL*power(rho,K,mod)*inv%mod!=m2)throw runtime_error("long entry phase");entries[key(z,m2,gap+2)].push_back({K+a,K+1});++longtests;}
   ++entrytests;m=1LL*m*rho%mod;
  }
 }
 vector<Node>nodes;nodes.reserve(100000);unordered_map<uint64_t,int>ids;ids.reserve(100000);auto add=[&](int z,int m,int b){uint64_t k=key(z,m,b);auto it=ids.find(k);if(it!=ids.end())return it->second;int i=nodes.size();if(i>=6000000)throw runtime_error("incomplete six-million node guard");ids[k]=i;nodes.push_back(Node{z,m,b,{},{},0});return i;};
 vector<int>roots;for(int b=0;b<=B;b++)roots.push_back(add(0,1,b));ll scanned=0,edges=0;
 for(size_t pos=0;pos<nodes.size();pos++){
  int z=nodes[pos].z,m=nodes[pos].m,b=nodes[pos].b,e=h[z],mod=p[max(R,e)],target=value(z,m),reward=min<int>(e,h[m-1]);vector<pair<int,int>>pred;
  for(int a=1;a<=64;a++){int want=target-8*reward+5*a;if(want-b>hi||want<lo)continue;ll zz=(1LL*p2[a]*(z-1)+2)%q;if(zz<0)zz+=q;int z0=zz*inv3%q,mod0=p[max(R,int(h[z0]))],common=min(mod,mod0);int residue=(ll(m)*p2[a]%common)*inv3%common;
   for(int m0=residue;m0<mod0;m0+=common){++scanned;int gap=want-value(z0,m0);if(gap>=0&&gap<=b)pred.push_back({add(z0,m0,b-gap),a});}
  }
  edges+=pred.size();nodes[pos].pred=move(pred);if(pos%100000==0)cerr<<"processed="<<pos<<" nodes="<<nodes.size()<<" scanned="<<scanned<<"\n";
 }
 vector<uint8_t>color(nodes.size(),0);ll processed=0;for(int root:roots){if(color[root])continue;vector<pair<int,size_t>>stack{{root,0}};color[root]=1;
  while(!stack.empty()){int n=stack.back().first;size_t &j=stack.back().second;if(j<nodes[n].pred.size()){int c=nodes[n].pred[j++].first;if(color[c]==1)throw runtime_error("augmented cycle");if(!color[c]){color[c]=1;stack.push_back({c,0});}continue;}
   Node &N=nodes[n];auto it=entries.find(key(N.z,N.m,N.b));if(it!=entries.end()){N.pairs=it->second;N.paths=it->second.size();}
   for(auto [c,a]:N.pred){if(color[c]!=2)throw runtime_error("topology");if(numeric_limits<uint64_t>::max()-N.paths<nodes[c].paths)throw runtime_error("path overflow");N.paths+=nodes[c].paths;for(auto [A,L]:nodes[c].pairs)N.pairs.push_back({A+a,L+1});}
   sort(N.pairs.begin(),N.pairs.end());N.pairs.erase(unique(N.pairs.begin(),N.pairs.end()),N.pairs.end());if(N.pairs.size()>1000000)throw runtime_error("pair guard");color[n]=2;++processed;stack.pop_back();
  }
 }
 if(processed!=(ll)nodes.size())throw runtime_error("unvisited graph");ll retainedentries=0;for(auto &N:nodes){auto it=entries.find(key(N.z,N.m,N.b));if(it!=entries.end())retainedentries+=it->second.size();}
 ll endpoint=1;for(int i=0;i<END;i++)endpoint*=5;
 cout<<"{\"status\":\"complete near-tight graph\",\"E\":"<<E<<",\"R\":"<<R<<",\"budget\":"<<B<<",\"endpoint_exponent\":"<<END<<",\"terminal\":"<<value(0,1)<<",\"minimum_potential\":"<<lo<<",\"maximum_potential\":"<<hi<<",\"nodes\":"<<nodes.size()<<",\"edges\":"<<edges<<",\"scanned_candidates\":"<<scanned<<",\"entries_checked\":"<<entrytests<<",\"long_entry_controls\":"<<longtests<<",\"retained_entries\":"<<retainedentries<<",\"acyclic\":true,\"records\":[";
 for(size_t i=0;i<roots.size();i++){if(i)cout<<",";Node &N=nodes[roots[i]];cout<<"{\"slack\":"<<N.b<<",\"score\":"<<-value(0,1)+N.b<<",\"abstract_paths\":"<<N.paths<<",\"pairs\":[";bool first=true;vector<Pair>survive;for(auto [A,L]:N.pairs){if(!(1<=L&&L<A))throw runtime_error("word clock");if(!first)cout<<",";first=false;cout<<"["<<A<<","<<L<<"]";if(power(3,L,endpoint)==power(2,A,endpoint))survive.push_back({A,L});}cout<<"],\"multiplier_survivors\":[";first=true;for(auto [A,L]:survive){if(!first)cout<<",";first=false;cout<<"["<<A<<","<<L<<"]";}cout<<"]}";}
 cout<<"]}\n";return 0;
 }catch(const exception&e){cerr<<e.what()<<"\n";return 1;}}
