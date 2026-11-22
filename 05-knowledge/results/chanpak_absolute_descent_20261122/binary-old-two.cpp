#include<bits/stdc++.h>
using namespace std;using I=long long;using R=vector<int>;
R flip(R a,int u,int v){a[u]^=1<<v;a[v]^=1<<u;return a;}
R extend(R a,int s){for(int i=0;i<(int)a.size();i++)if(!(s>>i&1))a[i]|=1<<a.size();a.push_back(s);return a;}
void dp(const R&a,array<array<I,10>,1024>& f){
 int n=a.size(),m=1<<n;
 for(int s=0;s<m;s++)f[s].fill(0);
 for(int i=0;i<n;i++)f[1<<i][i]=1;
 for(int s=1;s<m-1;s++)for(int i=0;i<n;i++)if(f[s][i])for(int t=a[i]&~s;t;t&=t-1){int j=__builtin_ctz(t);f[s|1<<j][j]+=f[s][i];}
}
I hp(const R&a){
 static array<array<I,10>,1024> f;
 dp(a,f); return accumulate(f[(1<<a.size())-1].begin(),f[(1<<a.size())-1].end(),0LL);
}
vector<I> all_apex(const R&a){
 int n=a.size(),m=1<<n;R b(n);
 for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(a[i]>>j&1)b[j]|=1<<i;
 static array<array<I,10>,1024> E,B;
 dp(a,E); dp(b,B);
 I w[10][10]{},L[10]{};
 for(int s=1;s<m-1;s++)for(int u=s;u;u&=u-1){
  int i=__builtin_ctz(u);if(!E[s][i])continue;
  for(int v=(m-1)^s;v;v&=v-1){int j=__builtin_ctz(v);w[i][j]+=E[s][i]*B[(m-1)^s][j];}
 }
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){w[i][j]+=w[j][i];w[j][i]=w[i][j];L[(a[i]>>j&1)?i:j]+=w[i][j];}
 vector<I> ans(m);
 ans[0]=accumulate(E[m-1].begin(),E[m-1].end(),0LL);
 for(int s=1;s<m;s++){
  int i=__builtin_ctz(s),q=s&(s-1);I v=ans[q]+L[i];
  for(int t=q;t;t&=t-1)v-=w[i][__builtin_ctz(t)];
  ans[s]=v;
 }
 return ans;
}

int main(){
 R P={22,44,88,49,98,69,11},B=flip(P,0,1);
 set<int> TARGET={377,379,753,755,757,759}; set<int> printed; I mixed_min=LLONG_MAX;
 for(int u=0;u<7;u++)for(int v=u+1;v<7;v++){
   R a=(u==0&&v==1)?B:flip(B,u,v);
   I mm=LLONG_MAX;
   for(int s=0;s<128;s++){
      R b=extend(a,s);auto vals=all_apex(b);
      for(int t=0;t<(int)vals.size();t++){
        if(0<s&&s<127&&0<(t&127)&&(t&127)<127) mm=min(mm,vals[t]);
        if(TARGET.contains(vals[t])&&!printed.contains(vals[t])){
           assert(hp(extend(b,t))==vals[t]);printed.insert(vals[t]);
           cout<<"EXACT "<<vals[t]<<" old flips (0,1) ";
           if(!(u==0&&v==1))cout<<"("<<u<<","<<v<<") ";
           cout<<"tags "<<s<<" "<<t<<endl;
        }
      }
   }if(u||v!=1)mixed_min=min(mixed_min,mm);cout<<"mixed-min "<<u<<" "<<v<<" "<<mm<<endl;
 }
 assert(!printed.contains(377)); assert(mixed_min==333);
}
