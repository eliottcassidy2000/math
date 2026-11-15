#include <bits/stdc++.h>
using namespace std;
struct Data{int n,F;vector<int> h,S,d,delta; Data(vector<int> r):n(r.size()),F((1<<n)-1),h(F+1),S(n),d(n),delta(n){
 int a[2048][11];memset(a,0,sizeof(a[0])*(F+1));h[0]=1;
 for(int m=1;m<=F;m++)for(int j=m;j;j&=j-1){
  int i=__builtin_ctz(unsigned(j)),q=m^(1<<i);
  if(!q)a[m][i]=1;
  else for(int u=q&r[i];u;u&=u-1)a[m][i]+=a[q][__builtin_ctz(unsigned(u))];
  h[m]+=a[m][i];
 }
 for(int m=1;m<=F;m++)for(int z=m;z;z&=z-1){
  int i=__builtin_ctz(unsigned(z)),x=a[m][i]*h[F^m];
  S[i]+=x;
  d[i]+=__builtin_popcount(unsigned(m))%2?x:-x;
 }
 for(int i=0;i<n;i++)delta[i]=S[i]-d[i]-h[F];
}};
int main(){
 vector<int> t={0},w;
 map<int,int> first;map<int,vector<int>> wit;
 int total=0;int literal=0;
 function<void(int)> rec=[&](int n){
  Data p(t); total++;
  int H=p.h.back();
  if(!first.count(H)||n<first[H]){first[H]=n;wit[H]=w;}
  if(n==10){
    for(int q=0;q<n;q++){
      total++;
      int k=H+p.S[q]-p.d[q];
      if(total%41999==0){
        auto r=t;
        for(int i=0;i<n;i++)if(i!=q)r[i]|=1<<n;
        r.push_back(1<<q);
        if(Data(r).h.back()!=k)throw runtime_error("literal leaf check");
        literal++;
      }
      if(k>=366|| (first.count(k)&&first[k]<=11)) continue;
      first[k]=11; wit[k]=w;wit[k].push_back(q);
    }
    return;
  }
  for(int q=0;q<n;q++){
    for(int i=0;i<n;i++)if(i!=q)t[i]|=1<<n;
    t.push_back(1<<q);w.push_back(q);rec(n+1);w.pop_back();t.pop_back();
    for(int i=0;i<n;i++)t[i]&=(1<<n)-1;
  }
 };
 rec(1);
 cout<<"checked "<<total<<" histories; "<<literal<<" independent literal final-layer DPs"<<endl;
 cout<<"holes below 366";for(int x=1;x<366;x+=2)if(!first.count(x))cout<<" "<<x;cout<<endl;for(int k:{23,33,35}){cout<<k<<" first "<<(first.count(k)?first[k]:0)<<" parents";for(auto x:wit[k])cout<<" "<<x;cout<<endl;}
}
