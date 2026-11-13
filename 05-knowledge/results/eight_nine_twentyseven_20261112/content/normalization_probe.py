"""Bounded physical same-H source changes, not an isomorphism census.

Compile the self-contained finite driver temporarily, then independently
recount every improvement with the existing Python subset-cover routine.
"""
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
CPP = r'''
#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <vector>
using namespace std;
using ll=long long;
using G=vector<int>;
ll fact[11];
vector<ll> counts(const G&t){
 int n=t.size(),s=1<<n;
 vector<array<ll,10>> d(s); vector<ll>h(s);h[0]=1;
 for(int i=0;i<n;i++)d[1<<i][i]=1;
 for(int m=1;m<s;m++)for(int i=0;i<n;i++)if(d[m][i]){
  h[m]+=d[m][i];int q=t[i]&~m;
  while(q){int b=q&-q;q-=b;d[m|b][__builtin_ctz(b)]+=d[m][i];}
 }
 return h;
}
ll content(const G&t,const vector<ll>&h){
 int n=t.size(),s=1<<n; vector<array<ll,11>>p(s);p[0][0]=1;
 for(int m=1;m<s;m++){
  int bit=m&-m;
  for(int sub=m;sub;sub=(sub-1)&m)if(sub&bit){
   int rem=m^sub;
   for(int k=0;k<=__builtin_popcount((unsigned)rem);k++)
    p[m][k+1]+=h[sub]*p[rem][k];
  }
 }
 ll g=0;for(int k=1;k<=n;k++)g=gcd(g,p[s-1][k]*fact[k]*fact[k-1]);
 return g;
}
void arrayjson(const G&a){cout<<'[';for(int i=0;i<(int)a.size();i++){if(i)cout<<',';cout<<a[i];}cout<<']';}
G flip(G t,int i,int j){t[i]^=1<<j;t[j]^=1<<i;return t;}
G pad(const G&t,int a,int b){
 int n=t.size(),N=n+a+b;G u(N);
 for(int i=0;i<N;i++)for(int j=i+1;j<N;j++){
  bool forward=!(i>=a&&j<a+n)||((t[i-a]>>(j-a))&1);
  if(forward)u[i]|=1<<j;else u[j]|=1<<i;
 }
 return u;
}
int main(){
 fact[0]=1;for(int i=1;i<=10;i++)fact[i]=i*fact[i-1];
 G control={56,117,129,22,100,76,141,59},spine={1,0,4,5,6,2,7,3},reset(8);
 for(int i=0;i<8;i++)reset[i]=255^control[i]^(1<<i);
 for(int j=1;j<8;j++)reset=flip(reset,spine[j-1],spine[j]);
 if(counts(control).back()!=513||content(control,counts(control))!=27||
    counts(reset).back()!=513||content(reset,counts(reset))!=3)return 2;
 cout<<"{\"type\":\"positive_control\",\"H\":513,\"g\":3,\"rows\":";arrayjson(reset);
 cout<<",\"base\":";arrayjson(control);cout<<",\"args\":";arrayjson(spine);cout<<"}\n";
 vector<pair<ll,G>> seeds={{55,{54,124,56,113,96,64,5}},
 {171,{42,28,25,112,97,70,7}},{175,{14,28,56,112,97,67,7}},
 {189,{22,44,88,49,98,69,11}}};
 for(auto entry:seeds){
  ll H=entry.first,target=gcd(H,3LL),best=content(entry.second,counts(entry.second));
  ll tries=0,neutral=0;set<G>seen;array<ll,4>bykind{},neutralbykind{};map<ll,ll>ghist;
  bool done=false;
  auto test=[&](const G&t,const G&base,int kind,const G&args){
   if(done)return;tries++;bykind[kind]++;
   auto h=counts(t);if(h.back()!=H)return;neutral++;neutralbykind[kind]++;
   if(!seen.insert(t).second)return;
   ll g=content(t,h);ghist[g]++;if(g>=best)return;best=g;
   cout<<"{\"type\":\"improvement\",\"H\":"<<H<<",\"g\":"<<g
       <<",\"kind\":"<<kind<<",\"rows\":";arrayjson(t);
   cout<<",\"base\":";arrayjson(base);cout<<",\"args\":";arrayjson(args);cout<<"}\n"<<flush;
   if(best==target)done=true;
  };
  for(int n=7;n<=10&&!done;n++)for(int a=0;a<=n-7&&!done;a++){
   G base=pad(entry.second,a,n-7-a);vector<pair<int,int>> edges;
   for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)edges.push_back({i,j});
   for(auto e:edges)if(!done)test(flip(base,e.first,e.second),base,0,{e.first,e.second});
   for(int u=0;u<(int)edges.size()&&!done;u++)for(int v=u+1;v<(int)edges.size()&&!done;v++){
    auto e=edges[u],f=edges[v];test(flip(flip(base,e.first,e.second),f.first,f.second),base,1,{e.first,e.second,f.first,f.second});
   }
   G path;
   function<void(int)>dfs=[&](int used){
    if(done)return;
    if((int)path.size()==n){
     G t(n);for(int i=0;i<n;i++)t[i]=((1<<n)-1)^base[i]^(1<<i);
     for(int j=1;j<n;j++)t=flip(t,path[j-1],path[j]);
     test(t,base,2,path);return;
    }
    int q=path.empty()?((1<<n)-1):(base[path.back()]&~used);
    while(q&&!done){int bit=q&-q;q-=bit;path.push_back(__builtin_ctz(bit));dfs(used|bit);path.pop_back();}
   };dfs(0);
   if(n<=9)for(int mask=0;mask<(1<<n)&&!done;mask++){
    G t=base;t.push_back(mask);for(int i=0;i<n;i++)if(!(mask>>i&1))t[i]|=1<<n;
    test(t,base,3,{mask});
   }
  }
  cout<<"{\"type\":\"summary\",\"H\":"<<H<<",\"minimum_found\":"<<best
      <<",\"parity_target\":"<<target<<",\"attempts\":"<<tries<<",\"neutral_attempts\":"<<neutral
      <<",\"distinct_neutral_labelled_graphs\":"<<seen.size()<<",\"attempts_by_kind\":";
  arrayjson(G(bykind.begin(),bykind.end()));cout<<",\"neutral_by_kind\":";
  arrayjson(G(neutralbykind.begin(),neutralbykind.end()));cout<<",\"distinct_neutral_g_histogram\":{";
  bool first=true;for(auto item:ghist){if(!first)cout<<',';first=false;cout<<'"'<<item.first<<"\":"<<item.second;}
  cout<<"}}\n"<<flush;
 }
}
'''


def main():
    with tempfile.TemporaryDirectory(prefix='.normalization-', dir=HERE) as tmp:
        source=Path(tmp)/'probe.cpp'; binary=Path(tmp)/'probe'
        source.write_text(CPP)
        subprocess.run(['c++','-std=c++17','-O3',str(source),'-o',str(binary)],check=True)
        run=subprocess.run([str(binary)],capture_output=True,text=True,check=True,timeout=180)
    lines=[json.loads(line) for line in run.stdout.splitlines()]
    spec=importlib.util.spec_from_file_location('content_check',HERE/'check.py')
    check=importlib.util.module_from_spec(spec);spec.loader.exec_module(check)
    for item in lines:
        if item['type'] not in ('improvement','positive_control'):
            continue
        c=check.covers_dp(item['rows'])
        assert c[0]==item['H'] and check.content(c)['contextual']==item['g']
        item['ordered_covers']=c
        item['polynomial_coefficients']=check.from_covers(c)
        if item['type']=='positive_control':
            base=item['base'];path=item['args'];n=len(base)
            reset=[((1<<n)-1)^row^(1<<i) for i,row in enumerate(base)]
            assert all(base[a]>>b&1 for a,b in zip(path,path[1:]))
            for a,b in zip(path,path[1:]):
                reset[a]^=1<<b;reset[b]^=1<<a
            assert reset==item['rows']
    receipt=dict(status='FINITE-EXACT bounded operation search; not global normalization',
                 universe='Four named order7 sources; every transitive prefix/suffix padding through order10; one and two edge reversals, every Hamilton-spine complement, and all one-vertex orientations from bases of order<=9. Stop each source at parity minimum.',
                 kinds=['one_edge_reversal','two_edge_reversals','native_spine_complement','one_vertex_extension'],
                 records=lines)
    (HERE/'normalization_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps([x for x in lines if x['type']=='summary'],indent=2))


if __name__=='__main__':
    main()
