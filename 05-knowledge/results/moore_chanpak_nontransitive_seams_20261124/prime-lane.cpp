#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <vector>
using namespace std;
using ll=long long;using Rows=vector<int>;using P=array<ll,3>;
P prof(const Rows&r){
 int n=r.size(),U=(1<<n)-1;
 vector<array<int,13>> dp(U+1);vector<int> h(U+1),q(U+1);
 for(int i=0;i<n;i++)dp[1<<i][i]=1;
 for(int s=1;s<=U;s++)for(int i=s;i;i&=i-1){
  int j=__builtin_ctz(i),a=dp[s][j];h[s]+=a;
  if(a)for(int x=r[j]&~s;x;x&=x-1)dp[s|(x&-x)][__builtin_ctz(x)]+=a;
 }
 for(int s=1;s<=U;s++)for(int x=(s-1)&s;x;x=(x-1)&s)if(x<(s^x))q[s]+=2*h[x]*h[s^x];
 ll t=0;for(int s=1;s<U;s++)t+=(ll)q[s]*h[U^s];
 return {h[U],q[U],t};
}
Rows add(const Rows&r,int x){Rows a=r;int n=r.size();for(int i=0;i<n;i++)if(!(x>>i&1))a[i]|=1<<n;a.push_back(x);return a;}
ll height(const P&a,const P&b){static constexpr int C[3][3]={{109,117,32},{94,117,40},{24,36,16}};ll h=0;for(int i=0;i<3;i++)for(int j=0;j<3;j++)h+=C[i][j]*a[i]*b[j];return h;}
bool prime(ll v){if(v<=1)return false;for(ll j=2;j*j<=v;j+=(j==2?1:2))if(v%j==0)return false;return true;}
void print(const P&p,const Rows&r){for(ll x:p)cout<<x<<',';cout<<" / ";for(int x:r)cout<<x<<',';}
int main(int argc,char**argv){
 int ord=argc>1?atoi(argv[1]):9, rounds=argc>2?atoi(argv[2]):200000;
 ifstream f(argc>3?argv[3]:"moore_n8_all.data");if(!f)return 2;int n;vector<Rows> db[9];
 while(f>>n){Rows a(n);for(int&i:a)f>>i;db[n].push_back(a);}
 map<P,Rows> Ys,Xs;
 for(int n=2;n<=6;n+=2)for(const auto&r:db[n]){P p=prof(r);if(p[0]%3)Ys.try_emplace(p,r);}
 for(const auto&r:db[7])for(int x=0;x<128;x++){auto t=add(r,x);P p=prof(t);if(p[0]%3)Ys.try_emplace(p,t);}
 struct Entry{P p,c;Rows r;};vector<Entry> bank;
 for(auto &[p,r]:Ys){auto [d,e,f]=p; bank.push_back({p,{2*d-17*e-8*f,17*d-4*f,8*d+4*e},r});}
 cout<<"Y bank "<<bank.size()<<endl;
 ll checked=0,hits=0, mod3hits=0, actual=0;
 auto scan=[&](const Rows&r){
  auto a=prof(r);if(a[0]%3==0)return;
  if(!Xs.try_emplace(a,r).second)return;checked++;
  for(auto &b:bank){
   ll d=a[0]*b.c[0]+a[1]*b.c[1]+a[2]*b.c[2];
   if(d!=2&&d!=-2)continue;
   hits++;ll h=height(a,b.p);
   if(h%3==0||(h+d)%3==0)continue;
   mod3hits++;
   cout<<"UNIT "<<h<<" "<<h+d<<" : ";print(a,r);cout<<" vs ";print(b.p,b.r);cout<<endl;
   if(prime(h)&&prime(h+d)){
    cout<<"PRIME "<<h<<" "<<h+d<<" : ";print(a,r);cout<<" vs ";print(b.p,b.r);cout<<endl;actual++;
   }
  }
 };
 for(int n=1;n<=7;n+=2)for(auto&r:db[n])scan(r);
 if(ord==9) for(auto &r:db[8])for(int x=0;x<256;x++)scan(add(r,x));
 mt19937_64 RNG(101);
 for(int t=0;t<rounds;t++){
  Rows r;
  if(t%2==0)r=db[7][RNG()%db[7].size()]; // isomorphism-uniform seven-base, then two random apices
  for(int i=r.size();i<ord;i++)r=add(r,RNG()%(1<<i)); // odd fixed final size
  scan(r);
  if(t%100000==99999)cout<<"progress "<<t+1<<" "<<checked<<" "<<hits<<" "<<mod3hits<<" "<<actual<<endl;
 }
 if(ord==9&&rounds==0)assert(checked==32820&&bank.size()==1103&&hits==1&&mod3hits==1&&actual==0);
 cout<<"RESULT "<<ord<<" "<<rounds<<" "<<checked<<" "<<hits<<" "<<mod3hits<<" "<<actual<<endl;
}
