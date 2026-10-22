#include <algorithm>
#include <array>
#include <iostream>
#include <vector>
using namespace std; using Sig=array<long long,4>;
struct Type{int k; Sig s;}; vector<Type> a; vector<Sig> v;
void rec(int idx,int rem,Sig s){if(!rem){v.push_back(s);return;}for(int i=idx;i<(int)a.size();i++){auto t=a[i];if(t.k>rem)break; Sig z=s;for(int j=0;j<4;j++)z[j]+=t.s[j];rec(i,rem-t.k,z);}}
int main(){for(long long k=1;k<=43;k++){a.push_back({(int)k,{k*(k+1),-k*(k+1)*(2*k-1),4*k*(k+1)*(5*k*k-7*k+3),-2*k*(k+1)*(14*k*k*k-37*k*k+39*k-15)}});if(k>=4)a.push_back({(int)k,{2*k*(k-1),-2*k*(k-1)*(4*k-7),8*k*(k-1)*(20*k*k-74*k+69),-4*k*(k-1)*(112*k*k*k-672*k*k+1378*k-963)}});}
a.push_back({6,{72,-1512,152928,-2011824}});a.push_back({7,{126,-4158,669816,-14196924}});a.push_back({8,{240,-13680,3844800,-143445600}});sort(a.begin(),a.end(),[](Type x,Type y){return x.k<y.k;});long long total=0;for(int n=1;n<=43;n++){v.clear();rec(0,n,{0,0,0,0});total+=v.size();sort(v.begin(),v.end());int dup=0;for(size_t i=1;i<v.size();i++)if(v[i]==v[i-1]){dup++;cout<<"collision rank "<<n<<" signature ";for(auto x:v[i])cout<<x<<",";cout<<"\n";}cout<<n<<" "<<v.size()<<" duplicates="<<dup<<" cumulative="<<total<<endl;}}
