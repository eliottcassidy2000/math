// Exact port of zero_ready_component_returns_20261103/interface/gap.count.
// Domain n<=30, 1<=bins<=3; feedback deletion must be transitive.
// Every intermediate weight counts distinct partial assignments, injectively
// encoded by a permutation of p+m<=n vertices and k ordered possibly empty
// bins. Thus weight<=n!*C(n+k-1,k-1)<=30!*496<2^117. The accepted total has
// the same injection. Unsigned128 therefore cannot overflow, intermediates
// included. This is a domain guard, not a practical factorial-time promise.
#include <algorithm>
#include <cstdint>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;
using W=unsigned __int128;
using S=uint64_t;
bool force_sparse=false;
void need(bool ok,const char* message){if(!ok)throw runtime_error(message);}
string decimal(W x){if(!x)return "0";string s;while(x){s.push_back('0'+x%10);x/=10;}reverse(s.begin(),s.end());return s;}
struct Gap{int left,right;bool forbidden;};
W count(const vector<uint32_t>& rows,vector<int> ports,int bins,const vector<int>& starts,const vector<int>& ends){
 const int n=rows.size();
 need(n<=30&&bins>=1&&bins<=3,"n<=30 and 1<=bins<=3 required");
 for(int v=0;v<n;v++){
  need((rows[v]>>n)==0&&((rows[v]>>v)&1)==0,"row bounds or diagonal");
  for(int u=0;u<v;u++)need(((rows[v]>>u)&1)+((rows[u]>>v)&1)==1,"not a tournament");
 }
 for(int v:ports)need(v>=0&&v<n,"feedback vertex");
 for(int b=0;b<bins;b++)for(int v:{starts[b],ends[b]})if(v>=0){need(v<n,"endpoint vertex");ports.push_back(v);}
 sort(ports.begin(),ports.end());ports.erase(unique(ports.begin(),ports.end()),ports.end());
 vector<int> chain;
 for(int v=0;v<n;v++)if(!binary_search(ports.begin(),ports.end(),v))chain.push_back(v);
 vector<int> degree(n);
 for(int v:chain)for(int u:chain)degree[v]+=(rows[v]>>u)&1;
 sort(chain.begin(),chain.end(),[&](int a,int b){return degree[a]!=degree[b]?degree[a]>degree[b]:a<b;});
 for(size_t i=0;i<chain.size();i++)for(size_t j=i+1;j<chain.size();j++)need((rows[chain[i]]>>chain[j])&1,"feedback deletion not transitive");
 const int p=ports.size(),g=p+bins;
 vector<S> powers(g+1,1);for(int i=1;i<=g;i++)powers[i]=3*powers[i-1]; // g<=33; 3^33<2^64.
 W total=0;vector<int> sizes(bins);
 vector<W> current,next;vector<S> active,next_active;
 const bool dense=g<=12&&!force_sparse;
 if(dense){current.resize(powers[g]);next.resize(powers[g]);}
 auto evaluate=[&](){
  vector<Gap> gaps;int pos=0;
  for(int b=0;b<bins;b++){
   int z=sizes[b];
   if((starts[b]>=0&&(!z||ports[pos]!=starts[b]))||(ends[b]>=0&&(!z||ports[pos+z-1]!=ends[b])))return;
   for(int j=0;j<=z;j++)gaps.push_back({j?ports[pos+j-1]:-1,j<z?ports[pos+j]:-1,(j==0&&starts[b]>=0)||(j==z&&ends[b]>=0)});
   pos+=z;
  }
  auto valid=[&](S state){
   for(int j=0;j<g;j++){
    int digit=state/powers[j]%3;auto x=gaps[j];
    if(digit==2||(digit==0&&x.left>=0&&x.right>=0&&!((rows[x.left]>>x.right)&1)))return false;
   }return true;
  };
  if(dense){
   active={0};current[0]=1;
   for(int v:chain){
    next_active.clear();
    for(S state:active){W weight=current[state];
     for(int j=0;j<g;j++){auto x=gaps[j];if(x.forbidden)continue;int old=state/powers[j]%3;
      if(!old&&x.left>=0&&!((rows[x.left]>>v)&1))continue;
      int fresh=x.right<0||((rows[v]>>x.right)&1)?1:2;
      S target=state+S(fresh)*powers[j]-S(old)*powers[j];
      if(!next[target])next_active.push_back(target);next[target]+=weight;
     }
    }
    for(S state:active)current[state]=0;
    current.swap(next);active.swap(next_active);
   }
   for(S state:active){if(valid(state))total+=current[state];current[state]=0;}
  }else{
   unordered_map<S,W> states{{0,1}},following;
   for(int v:chain){following.clear();
    for(auto [state,weight]:states)for(int j=0;j<g;j++){auto x=gaps[j];if(x.forbidden)continue;int old=state/powers[j]%3;
     if(!old&&x.left>=0&&!((rows[x.left]>>v)&1))continue;
     int fresh=x.right<0||((rows[v]>>x.right)&1)?1:2;
     following[state+S(fresh)*powers[j]-S(old)*powers[j]]+=weight;
    }states.swap(following);
   }for(auto [state,weight]:states)if(valid(state))total+=weight;
  }
 };
 function<void(int,int)> compose=[&](int b,int left){if(b==bins-1){sizes[b]=left;evaluate();return;}for(int z=0;z<=left;z++){sizes[b]=z;compose(b+1,left-z);}};
 do{compose(0,p);}while(next_permutation(ports.begin(),ports.end()));
 return total;
}
int main(int argc,char** argv){ios::sync_with_stdio(false);cin.tie(nullptr);
 try{need(argc==1||(argc==2&&string(argv[1])=="--sparse"),"only --sparse option is supported");force_sparse=argc==2;
 int n,k,nf,ns,ne;while(cin>>n>>k>>nf>>ns>>ne){
  need(n>=0&&n<=30&&k>=1&&k<=3&&nf>=0&&nf<=30&&ns>=0&&ns<=k&&ne>=0&&ne<=k,"invalid header");
  vector<uint32_t> rows(n);for(auto& x:rows)need(bool(cin>>x),"missing rows");
  vector<int> feedback(nf);for(auto&x:feedback)need(bool(cin>>x),"missing feedback");
  vector<int> starts(k,-1),ends(k,-1);
  for(int kind=0;kind<2;kind++)for(int i=0;i<(kind?ne:ns);i++){int b,v;need(bool(cin>>b>>v),"missing endpoint");need(b>=0&&b<k&&v>=0&&v<n,"endpoint bounds");auto& a=kind?ends:starts;need(a[b]<0,"duplicate endpoint bin");a[b]=v;}
  cout<<decimal(count(rows,feedback,k,starts,ends))<<'\n'<<flush;
 }}catch(const exception& e){cerr<<e.what()<<'\n';return 2;}
}
