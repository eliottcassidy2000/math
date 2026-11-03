// Exact bounded search for substitution-compatible two-complement returns.
// No output is an existence proof without replaying the emitted literal program.
#include <algorithm>
#include <bit>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <vector>
using namespace std;
using I=int64_t; using T=vector<int>;
struct DP { vector<I> h; vector<vector<I>> e; };
DP dp(const T& t) {
 int n=t.size(),f=(1<<n)-1; DP z{vector<I>(f+1),vector<vector<I>>(f+1,vector<I>(n))};z.h[0]=1;
 for(int s=1;s<=f;s++)for(int i=0;i<n;i++)if(s>>i&1){
  int r=s^(1<<i);I v=!r;
  for(int j=0;j<n;j++)if((r>>j&1)&&(t[j]>>i&1))v+=z.e[r][j];
  z.e[s][i]=v;z.h[s]+=v;
 }return z;
}
T opposite(const T&t){int f=(1<<t.size())-1;T u=t;for(int i=0;i<(int)t.size();i++)u[i]=f^t[i]^(1<<i);return u;}
vector<T> paths(const T&t){
 auto s=dp(opposite(t)).e;int n=t.size();vector<T> out;T p;
 function<void(int)> rec=[&](int left){if(!left){out.push_back(p);return;}for(int v=0;v<n;v++)
  if((left>>v&1)&&s[left][v]&&(p.empty()||(t[p.back()]>>v&1))){p.push_back(v);rec(left^(1<<v));p.pop_back();}};
 rec((1<<n)-1);return out;
}
T native(const T&t,const T&p){T u=opposite(t);for(int i=1;i<(int)p.size();i++){u[p[i-1]]^=1<<p[i];u[p[i]]^=1<<p[i-1];}return u;}
T birth(T t,int mask){int n=t.size();for(int i=0;i<n;i++)if(mask>>i&1)t[i]|=1<<n;t.push_back(((1<<n)-1)^mask);return t;}
bool strong(const T&t){for(const T&u:{t,opposite(t)}){int reached=1,old=0;while(old!=reached){old=reached;for(int i=0;i<(int)t.size();i++)if(reached>>i&1)reached|=u[i];}if(reached!=(1<<(int)t.size())-1)return false;}return true;}
T clone(const T&t,int q){int n=t.size();T ad;for(int i=0;i<n;i++){ad.push_back(i);if(i==q)ad.push_back(i);}T out(n+1);
 for(int i=0;i<=n;i++)for(int j=0;j<=n;j++)if((ad[i]==ad[j]&&i<j)||(ad[i]!=ad[j]&&(t[ad[i]]>>ad[j]&1)))out[i]|=1<<j;return out;}
vector<I> ears(const T&t){
 int n=t.size(),f=(1<<n)-1;auto a=dp(t),b=dp(opposite(t));vector<vector<I>> w(n,vector<I>(n));
 for(int s=1;s<f;s++)for(int i=0;i<n;i++)if(a.e[s][i])for(int j=0;j<n;j++)if(b.e[f^s][j])w[i][j]+=a.e[s][i]*b.e[f^s][j];
 vector<I> l(n),out(f+1,a.h[f]);for(int i=0;i<n;i++){l[i]=a.e[f][i]-b.e[f][i];for(int j=0;j<n;j++)l[i]+=w[i][j];}
 for(int s=1;s<=f;s++){int i=countr_zero((unsigned)s),r=s^(1<<i);out[s]=out[r]+l[i];for(int j=0;j<n;j++)if(r>>j&1)out[s]-=w[i][j]+w[j][i];}return out;
}
void print(const T&t){cout<<'[';for(int i=0;i<(int)t.size();i++){if(i)cout<<',';cout<<t[i];}cout<<']';}
int main(int argc,char**argv){
 if(argc<3)return 2;ifstream in(argv[1]);int n;in>>n;T t(n);for(auto&v:t)in>>v;
 int limit=stoi(argv[2]),offset=argc>3?stoi(argv[3]):0,start=argc>4?stoi(argv[4]):0;I h=dp(t).h.back();vector<int> ports;vector<I> h2(n);
 for(int v=0;v<n;v++)if(popcount((unsigned)t[v])==1||popcount((unsigned)t[v])==n-2){ports.push_back(v);h2[v]=dp(clone(t,v)).h.back();}
 struct Candidate{I h;T p,t;};vector<Candidate> first;auto pp=paths(t);
 // An independent complete census in the zero/one-complement, one-birth language.
 vector<int> one(2),one_strong(2);vector<T> weak_targets;auto check1=[&](const T&u){auto a=ears(u);for(int m=0;m<(int)a.size();m++)for(int i=0;i<2;i++)if(a[m]==3*h+2+2*i){one[i]++;T w=birth(u,m);if(strong(w))one_strong[i]++;else weak_targets.push_back(w);}};
 check1(t);
 for(auto&p:pp){T u=native(t,p);I a=dp(u).h.back();if(a<=3*h+4)check1(u);first.push_back({a,p,u});}
 cerr<<"one-complement census "<<h<<" paths "<<pp.size()<<" all hits "<<one[0]<<' '<<one[1]<<" strong hits "<<one_strong[0]<<' '<<one_strong[1]<<'\n';
 if(!limit){cout<<"{\"source_H\":"<<h<<",\"paths\":"<<pp.size()<<",\"all_hits\":";print(one);cout<<",\"strong_hits\":";print(one_strong);cout<<",\"weak_targets\":[";for(int j=0;j<(int)weak_targets.size();j++){if(j)cout<<',';print(weak_targets[j]);}cout<<"]}\n";return 0;}
 sort(first.begin(),first.end(),[](auto&a,auto&b){return a.h!=b.h?a.h<b.h:a.p<b.p;});
 I examined=0,hits=0;
 for(int rank=start;rank<min(limit,(int)first.size());rank++){
  auto&a=first[rank];vector<int> prev(n,-1),next(n,-1);for(int i=1;i<n;i++){prev[a.p[i]]=a.p[i-1];next[a.p[i-1]]=a.p[i];}
  cerr<<"outer "<<rank<<" H "<<a.h<<'\n';
  for(auto&q:paths(a.t)){
   vector<int> pv(n,-1),nx(n,-1),eligible;for(int i=1;i<n;i++){pv[q[i]]=q[i-1];nx[q[i-1]]=q[i];}
   for(int v:ports)if(prev[v]==pv[v]&&next[v]==nx[v])eligible.push_back(v);
   if(eligible.empty())continue;T u=native(a.t,q);if(dp(u).h.back()>3*h+4)continue;
   auto counts=ears(u);examined++;
   for(int m=0;m<(1<<n);m++)if((counts[m]==3*h+2||counts[m]==3*h+4)&&(!offset||counts[m]==3*h+offset)){
    hits++;for(int v:eligible){
     bool indegree_one=popcount((unsigned)t[v])==n-2;
     if(bool(m>>v&1)!=indegree_one)continue;
     T target=birth(u,m);if(!strong(target))continue;I z=dp(clone(target,v)).h.back();
     if(z==3*h2[v]+counts[m]-3*h){
      cout<<"{\"source_H\":"<<h<<",\"port\":"<<v<<",\"offset\":"<<counts[m]-3*h<<",\"clone_source_H\":"<<h2[v]<<",\"clone_target_H\":"<<z<<",\"rank\":"<<rank<<",\"examined\":"<<examined<<",\"first_path\":";print(a.p);cout<<",\"second_path\":";print(q);cout<<",\"mask\":"<<m<<",\"target\":";print(target);cout<<"}\n";cout.flush();return 0;
     }
    }
   }
  }
 }
 cerr<<"bounded end examined "<<examined<<" scalar hits "<<hits<<'\n';return 0;
}
