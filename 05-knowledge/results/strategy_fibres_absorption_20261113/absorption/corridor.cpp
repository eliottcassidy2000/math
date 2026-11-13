// Exact labelled tournament census. No external graph data or assertions.
// Encoding: pairs (0,1),(0,2),(1,2),(0,3),...; bit 1 orients i -> j.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>
using namespace std;
struct DSU {
  vector<int> p, sz;
  explicit DSU(int n):p(n),sz(n,1){iota(p.begin(),p.end(),0);}
  int root(int x){while(p[x]!=x){p[x]=p[p[x]];x=p[x];}return x;}
  void join(int a,int b){a=root(a);b=root(b);if(a==b)return;if(sz[a]<sz[b])swap(a,b);p[b]=a;sz[a]+=sz[b];}
};
int hp(int n, unsigned code){
  unsigned out[8]={}; int k=0;
  for(int j=1;j<n;++j)for(int i=0;i<j;++i,++k)
    if(code>>k&1)out[i]|=1u<<j;else out[j]|=1u<<i;
  int dp[256][8]={};for(int i=0;i<n;++i)dp[1<<i][i]=1;
  for(int mask=1;mask<(1<<n);++mask)for(int v=0;v<n;++v)if(dp[mask][v]){
    unsigned todo=out[v]&~unsigned(mask);
    while(todo){int w=__builtin_ctz(todo);todo&=todo-1;dp[mask|(1<<w)][w]+=dp[mask][v];}
  }
  int ans=0;for(int i=0;i<n;++i)ans+=dp[(1<<n)-1][i];return ans;
}
int brute(int n,unsigned code){
  bool adj[8][8]={};int k=0;for(int j=1;j<n;++j)for(int i=0;i<j;++i,++k){adj[i][j]=code>>k&1;adj[j][i]=!adj[i][j];}
  vector<int>p(n);iota(p.begin(),p.end(),0);int ans=0;
  do{bool ok=true;for(int i=1;i<n;++i)ok&=adj[p[i-1]][p[i]];ans+=ok;}while(next_permutation(p.begin(),p.end()));return ans;
}
struct Range {int lo=1000000000,hi=0;long long size=0;set<int> values;};
int main(int argc,char**argv){
  int cap=argc>1?stoi(argv[1]):7;if(cap<1||cap>7)throw runtime_error("cap must lie in 1..7");
  int jump=argc>2?stoi(argv[2]):1;if(jump<1)throw runtime_error("positive normalized jump required");
  int neutral_radius=argc>3?stoi(argv[3]):1;if(neutral_radius<1||neutral_radius>3)throw runtime_error("neutral radius must lie in 1..3");
  bool pair_steps=argc>4&&stoi(argv[4])!=0;
  vector<int>offset(cap+2,0),count(cap+1);for(int n=1;n<=cap;++n){count[n]=1<<(n*(n-1)/2);offset[n+1]=offset[n]+count[n];}
  vector<int> H(offset[cap+1]);DSU dsu(H.size());
  long long checks=0, unit_edges=0, padding_edges=0;
  cout<<"{\n  \"encoding\": \"pairs increasing upper endpoint, bit1 lower to upper\",\n  \"normalized_jump_bound\":"<<jump<<",\n  \"neutral_arc_radius\":"<<neutral_radius<<",\n  \"paired_step_normalized_jump\":"<<(pair_steps?jump:0)<<",\n  \"isomorphisms_identified\":true,\n  \"orders\": [\n";
  for(int n=1;n<=cap;++n){
    int m=n*(n-1)/2;set<int>spectrum;DSU iso(count[n]);
    for(int c=0;c<count[n];++c){int h=hp(n,c);if(h%2!=1)throw runtime_error("parity");H[offset[n]+c]=h;spectrum.insert(h);
      if(n<=5 || c%8191==0){if(h!=brute(n,c))throw runtime_error("independent permutation count");checks++;}
    }
    long long ne=0;for(int c=0;c<count[n];++c)for(int e=0;e<m;++e)if(!(c>>e&1)){
      int a=offset[n]+c,b=offset[n]+(c^(1<<e));if(abs(H[a]-H[b])<=2*jump){dsu.join(a,b);ne++;}
    }
    // Adjacent transpositions generate all relabellings. This also permits
    // deleting any source/sink, not just the last-labelled vertex.
    auto edge=[](int a,int b){if(a>b)swap(a,b);return b*(b-1)/2+a;};
    for(int a=0;a+1<n;++a){int b=a+1;vector<pair<int,int>>pairs;
      for(int v=0;v<n;++v)if(v!=a&&v!=b)pairs.emplace_back(edge(a,v),edge(b,v));
      for(int c=0;c<count[n];++c){int d=c^(1<<edge(a,b));
        for(auto[p,q]:pairs){int delta=((c>>p)^(c>>q))&1;d^=(delta<<p)|(delta<<q);}
        if(H[offset[n]+c]!=H[offset[n]+d])throw runtime_error("isomorphism invariant");
        dsu.join(offset[n]+c,offset[n]+d);
        iso.join(c,d);
      }
    }
    unit_edges+=ne;
    long long neutral_pair_edges=0;
    if(neutral_radius>=2)for(int c=0;c<count[n];++c)for(int e=0;e<m;++e)for(int f=e+1;f<m;++f){
      int d=c^(1<<e)^(1<<f);if(c<d&&abs(H[offset[n]+c]-H[offset[n]+d])<=(pair_steps?2*jump:0)){dsu.join(offset[n]+c,offset[n]+d);neutral_pair_edges++;}
    }
    long long triple_edges=0;vector<int>bridge5961;
    if(neutral_radius>=3)for(int c=0;c<count[n];++c)if(iso.root(c)==c)for(int e=0;e<m;++e)for(int f=e+1;f<m;++f)for(int g=f+1;g<m;++g){
      int d=c^(1<<e)^(1<<f)^(1<<g);
      if(abs(H[offset[n]+c]-H[offset[n]+d])<=(pair_steps?2*jump:0)){
        dsu.join(offset[n]+c,offset[n]+d);triple_edges++;
        if(bridge5961.empty()&&min(H[offset[n]+c],H[offset[n]+d])==59&&max(H[offset[n]+c],H[offset[n]+d])==61)bridge5961={c,d,e,f,g};
      }
    }
    if(n>1){int oldm=(n-1)*(n-2)/2;for(int c=0;c<count[n-1];++c)for(int kind=0;kind<2;++kind){
      int large=c|(kind?(((1<<(n-1))-1)<<oldm):0),a=offset[n-1]+c,b=offset[n]+large;
      if(H[a]!=H[b])throw runtime_error("neutral padding count");dsu.join(a,b);padding_edges++;
    }}
    map<int,Range> comps;for(int v=0;v<offset[n+1];++v){auto&r=comps[dsu.root(v)];r.lo=min(r.lo,H[v]);r.hi=max(r.hi,H[v]);r.size++;r.values.insert(H[v]);}
    vector<Range> ranges;for(auto &[id,r]:comps){if(jump==1&&int(r.values.size())!=(r.hi-r.lo)/2+1)throw runtime_error("connected unit component interval");ranges.push_back(r);}
    sort(ranges.begin(),ranges.end(),[](const Range&a,const Range&b){return a.lo<b.lo||(a.lo==b.lo&&a.hi>b.hi);});
    vector<pair<int,int>> maximal;for(auto&r:ranges){bool contained=false;for(auto q:maximal)if(q.first<=r.lo&&r.hi<=q.second){contained=true;break;}if(!contained)maximal.emplace_back(r.lo,r.hi);}
    int anchor_hi=0, anchor_components=0,anchor_gap=0;long long anchor_vertices=0;for(auto&r:ranges)if(r.values.count(23)){anchor_hi=max(anchor_hi,r.hi);anchor_components++;anchor_vertices+=r.size;int prev=*r.values.begin();for(int h:r.values){anchor_gap=max(anchor_gap,(h-prev)/2);prev=h;}}
    cout<<"    {\"n\":"<<n<<",\"labelled\":"<<count[n]<<",\"spectrum\":[";bool first=true;for(int h:spectrum){if(!first)cout<<",";first=false;cout<<h;}
    int isototal=0;for(int c=0;c<count[n];++c)isototal+=iso.root(c)==c;
    cout<<"],\"isomorphism_classes\":"<<isototal<<",\"bounded_flip_edges\":"<<ne<<",\"neutral_pair_edges\":"<<neutral_pair_edges<<",\"triple_representative_edges\":"<<triple_edges<<",\"bridge5961\":[";for(int i=0;i<int(bridge5961.size());++i){if(i)cout<<",";cout<<bridge5961[i];}
    cout<<"],\"cumulative_components\":"<<comps.size()<<",\"anchor23_max\":"<<anchor_hi<<",\"anchor23_components\":"<<anchor_components<<",\"anchor23_vertices\":"<<anchor_vertices<<",\"anchor23_max_gap\":"<<anchor_gap<<",\"maximal_component_hulls\":[";
    first=true;for(auto q:maximal){if(!first)cout<<",";first=false;cout<<"["<<q.first<<","<<q.second<<"]";}cout<<"]}"<<(n==cap?"\n":",\n");
    cerr<<"n="<<n<<" max="<<*spectrum.rbegin()<<" comps="<<comps.size()<<" anchor23_max="<<anchor_hi<<"\n";
  }
  cout<<"  ],\n  \"neutral_bridge25\":";
  bool found=false;
  if(jump==1&&neutral_radius==1&&cap>=7){
    map<int,pair<int,int>>bounds;for(int v=0;v<int(H.size());++v){int r=dsu.root(v);auto it=bounds.find(r);if(it==bounds.end())bounds[r]={H[v],H[v]};else{it->second.first=min(it->second.first,H[v]);it->second.second=max(it->second.second,H[v]);}}
    int n=7,m=21;for(int c=0;c<count[n]&&!found;++c)if(H[offset[n]+c]==25&&bounds[dsu.root(offset[n]+c)]==make_pair(23,25))
      for(int e=0;e<m&&!found;++e)for(int f=e+1;f<m&&!found;++f){int d=c^(1<<e)^(1<<f);
        if(H[offset[n]+d]==25&&bounds[dsu.root(offset[n]+d)]==make_pair(25,31)){
          cout<<"{\"n\":7,\"codes\":["<<c<<","<<d<<"],\"flipped_bits\":["<<e<<","<<f<<"],\"unit_component_intervals\":[[23,25],[25,31]]}";found=true;
        }
      }
  }
  if(!found)cout<<"null";
  cout<<",\n  \"independent_permutation_checks\":"<<checks<<",\n  \"bounded_flip_edges\":"<<unit_edges<<",\n  \"neutral_padding_edges\":"<<padding_edges<<",\n  \"status\":\"FINITE-EXACT\"\n}\n";
}
