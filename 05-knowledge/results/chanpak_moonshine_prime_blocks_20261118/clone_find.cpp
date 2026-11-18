#include <bits/stdc++.h>
using namespace std;
using ll=long long; using V=vector<int>;
pair<ll,ll> profile(const V& a) {
 int n=a.size(),f=(1<<n)-1;
 ll d[1<<14][14],h[1<<14]={}; for(int i=0;i<=f;i++)fill(d[i],d[i]+n,0);h[0]=1;
 for(int i=0;i<n;i++) d[1<<i][i]=1;
 for(int s=1;s<=f;s++) for(int x=s;x;x&=x-1){
   int i=__builtin_ctz(x); ll z=d[s][i]; if(!z)continue;h[s]+=z;
   for(int y=a[i]&~s;y;y&=y-1)d[s|(y&-y)][__builtin_ctz(y)]+=z;
 }
 ll b=0;for(int s=0;s<=f;s++)b+=h[s]*h[f^s];
 return {h[f],b};
}
V clone(const V &a,int k) {int x=a.size()-1, B=((1<<k)-1)<<x;V r;for(int i=0;i<x;i++)r.push_back((a[i]&~(1<<x))|((a[i]>>x&1)?B:0)); for(int i=0;i<k;i++)r.push_back(a.back()|(B&-(1<<(x+i+1))));return r;}
V cut(const V &a,int b) {auto r=a;int f=(1<<a.size())-1;for(int i=0;i<(int)a.size();i++)r[i]^=(b>>i&1)?f^b:b;return r;}
int main(int arg,char**argv){
 int n=atoi(argv[1]), z=atoi(argv[2]);bool input=arg>3;mt19937 random(103139);
 V base(n);int f=(1<<(n-1))-1,total=0,candidates=0;
 while(z--){fill(base.begin(),base.end(),0);if(input){for(int i=0;i<n-1;i++)if(!(cin>>base[i]))return 0;}
 else for(int i=0;i<n-1;i++)for(int j=i+1;j<n-1;j++)if(random()&1)base[i]|=1<<j;else base[j]|=1<<i;
 base.back()=f;total++;map<ll,vector<pair<int,ll>>> groups;
 for(int i=0;i<=f;i++){auto [h,b]=profile(cut(base,i));groups[b].push_back({i,h});}
 for(auto &[height,list]:groups)if(list.size()>1){
   map<vector<ll>, vector<pair<int,vector<ll>>>> pairgroups;
   for(auto [mask,h]:list){
     V A=cut(base,mask); auto [h2,b2]=profile(clone(A,2));
     pairgroups[{height,b2}].push_back({mask,{h,h2}});
   }
   for(auto &[key,pairs]: pairgroups) if(pairs.size()>1) {
     map<vector<ll>,pair<int,vector<ll>>> exact;
     for (auto [mask,H]:pairs){
       V A=cut(base,mask);
       int maxwidth=2+(n-1)/2;
       vector<ll> B=key;
       // Smaller cutoff r: interpolate exponential basis from width 1..r.
       // For these sizes fixed cutoff is simpler and still exact.
       for (int k=3;k<=maxwidth;k++){auto [h,b]=profile(clone(A,k));H.push_back(h);B.push_back(b);}
       auto [it,ok]=exact.emplace(B,make_pair(mask,H));
       if(ok) continue;
       candidates++;
       if(it->second.second!=H){
         cout<<"hit at input "<<total<<" masks "<<it->second.first<<" "<<mask<<" width "<<maxwidth<<"\n";
         for(auto t:base)cout<<t<<" ";cout<<endl;
         for(auto t:B)cout<<t<<" ";cout<<endl;
         for(auto t:it->second.second)cout<<t<<" ";cout<<endl;
         for(auto t:H)cout<<t<<" ";cout<<endl;
         return 0;
       }
     }
   }
 }
 }
 cout<<"no hit with exact buckets "<<total<<" checked pairs "<<candidates<<endl;
}
