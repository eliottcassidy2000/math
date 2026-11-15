#include <bits/stdc++.h>
using namespace std; using ul=uint64_t; using code=__uint128_t;
code canon(const array<int,16>& a,int n) {
  array<pair<ul,int>,16> v;
  for(int i=0;i<n;i++){
    int d=__builtin_popcount(a[i]);ul z=ul(d)<<60;
    for(int m=a[i];m;m&=m-1) z+=1ULL<<(4*__builtin_popcount(a[__builtin_ctz(m)]));
    v[i]={z,i};
  }
  sort(v.begin(),v.begin()+n);
  code best=~code(0);
  function<void(int)> calc=[&](int l) {
    if(l==n) {code x=0;int p=0;
      for(int i=0;i<n;i++) for(int j=i+1;j<n;j++) x|=code((a[v[i].second]>>v[j].second)&1)<<p++;
      best=min(best,x);return;
    }
    int r=l+1;while(r<n && v[r].first==v[l].first)++r;
    do {calc(r);} while(next_permutation(v.begin()+l,v.begin()+r));
  };
  calc(0);return best;
}
array<int,16> decode(code x,int n){
  array<int,16>a={};int p=0;
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)if((x>>p++)&1)a[i]|=1<<j;else a[j]|=1<<i;
  return a;
}
pair<ul,ul> profile(const array<int,16>& a,int n){
  int full=(1<<n)-1;vector<array<ul,15>> dp(1<<n); vector<ul> h(1<<n); h[0]=1;
  for(int i=0;i<n;i++)dp[1<<i][i]=1;
  for(int s=1;s<=full;s++)for(int w=s;w;w&=w-1){
    int i=__builtin_ctz(w); ul t=dp[s][i]; h[s]+=t;
    if(t)for(int m=a[i]&~s;m;m&=m-1){int j=__builtin_ctz(m);dp[s|1<<j][j]+=t;}
  }
  ul p=0;for(int s=1;s<full;s++)p+=h[s]*h[full^s];
  return {h[full],p};
}
bool strong(const array<int,16>& a,int n){
 int full=(1<<n)-1;
 for(int v=0;v<n;v++){
   int reached=1<<v,pending=reached;
   while(pending){int u=__builtin_ctz(pending);pending&=pending-1;
   int fresh=a[u]&~reached;pending|=fresh;reached|=fresh;}
   if(reached!=full)return false;
 }
 return true;
}
string code_str(code x) {
    if(x==0) return "0";
    string s;
    while(x){s+=char('0'+x%10);x/=10;}
    reverse(s.begin(),s.end());
    return s;
}
struct chash {
   size_t operator()(code x) const {return (ul)x ^ rot((ul)(x>>64));}
   static ul rot(ul x){return (x<<29)|(x>>35);}
};
using M=unordered_map<code,pair<ul,ul>,chash>;
int main(int argc,char**argv) {
 if(argc!=3){cerr<<"usage: capped FINAL_ORDER Q_CAP\n";return 2;}
 int target=stoi(argv[1]), cap=stoi(argv[2]);
 if(target<1||target>15||cap<(1<<target)) return 2;
 M old={{0,{1,0}}};
 for(int n=2;n<=target;++n) {
    M next;
    ul qcap=cap>>(target-n);
    vector<pair<code,pair<ul,ul>>> input(old.begin(),old.end());
    #pragma omp parallel
    {
      M local;
      #pragma omp for schedule(dynamic,128)
      for(size_t k=0;k<input.size();k++) {
        auto [parent,hp]=input[k]; auto [h,p]=hp;
        auto a=decode(parent,n-1);
        int full=(1<<(n-1))-1;
        for(int m=0;m<=full;m++) {
           int d=__builtin_popcount(unsigned(m));
           bool skip=false;
           for(int i=0;i<n-1;i++) if(__builtin_popcount(unsigned(a[i]))+!(m>>i&1)<d){skip=true;break;}
           if(skip)continue;
           int tris=0;
           for(int w=m;w;w&=w-1) tris+=__builtin_popcount(unsigned(a[__builtin_ctz(unsigned(w))]&(full^m)));
           if(2*(2*h+p)+(ul(tris)<<(n-1))>qcap)continue;
           auto b=a; b[n-1]=m;
           for(int i=0;i<n-1;i++)if(!(m>>i&1)) b[i]|=1<<(n-1);
           code z=canon(b,n);
           if(local.count(z))continue;
           auto rr=profile(b,n);
           if(2*rr.first+rr.second<=qcap)local.emplace(z,rr);
        }
      }
      #pragma omp critical
      next.merge(local);
    }
    old.swap(next);
    cerr<<n<<" "<<qcap<<" "<<old.size()<<endl;
 }
 map<pair<ul,ul>,array<code,2>> witnesses;
 for(auto [z,hp]:old) {
    auto it=witnesses.try_emplace(hp,array<code,2>{~code(0),~code(0)}).first;
    it->second[0]=min(it->second[0],z);
    if(strong(decode(z,target),target))it->second[1]=min(it->second[1],z);
 }
 for(auto [hp,z]:witnesses)
   cout<<target<<" "<<hp.first<<" "<<hp.second<<" "<<code_str(z[0])<<" "<<(z[1]==~code(0)?"-":code_str(z[1]))<<"\n";
}
