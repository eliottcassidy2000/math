#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>
using namespace std;using ll=long long;
ll factor(ll x){if(x<2)return 0;for(ll d=2;d*d<=x;d++)if(x%d==0)return d;return 1;}
void need(bool x,const char*s){if(!x)throw runtime_error(s);}
int main(){map<ll,vector<int>> rows;map<ll,array<ll,3>> profiles;ll tested=0;vector<ll>accepted_orders;const ll LIMIT=73529109;
for(int n=1;n<=7;n++){int e=n*(n-1)/2;ll accepted=0;for(int code=0;code<(1<<e);code++){vector<int>a(n);int bit=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,bit++)a[(code>>bit&1)?i:j]|=1<<((code>>bit&1)?j:i);bool sorted=true;for(int i=1;i<n;i++)if(__builtin_popcount((unsigned)a[i-1])>__builtin_popcount((unsigned)a[i])){sorted=false;break;}if(!sorted)continue;accepted++;tested++;int M=1<<n;ll dp[128][7]={},h[128]={},c2[128]={};h[0]=1;for(int i=0;i<n;i++)dp[1<<i][i]=1;for(int S=1;S<M;S++){for(int v=0;v<n;v++){h[S]+=dp[S][v];for(int todo=a[v]&~S;todo;todo&=todo-1){int b=todo&-todo;dp[S|b][__builtin_ctz((unsigned)b)]+=dp[S][v];}}for(int U=(S-1)&S;U;U=(U-1)&S)c2[S]+=h[U]*h[S^U];c2[S]/=2;}int all=M-1;ll c3=0;for(int U=all;U;U=(U-1)&all)if(U&1)c3+=h[U]*c2[all^U];ll x=h[all],y=c2[all],z=c3;ll H=5*x*x+12*x*y+28*y*y+48*y*z+72*z*z;if(H<LIMIT&&!rows.count(H)){rows[H]=a;profiles[H]={x,y,z};}}
accepted_orders.push_back(accepted);}
ll np=0,n7=0;vector<ll>sevens;for(auto &[H,a]:rows){if(factor(H)==1)np++;if(H%7==0&&factor(H/7)==1){n7++;sevens.push_back(H);}}
need(accepted_orders==vector<ll>({1,1,3,9,61,529,11559}),"complete degree-filter counts");need(rows.size()==119&&np==8&&n7==3,"complete gate spectrum counts");need(sevens==vector<ll>({37050461,46937093,56228501}),"exact seven-prime list");
cout<<"{\n  \"status\": \"PASS\",\n  \"exclusive_limit\": "<<LIMIT<<",\n  \"maximum_input_order\": 7,\n  \"degree_sorted_counts\": [";for(int i=0;i<7;i++){if(i)cout<<",";cout<<accepted_orders[i];}cout<<"],\n  \"tested_inputs\": "<<tested<<",\n  \"distinct_values\": "<<rows.size()<<",\n  \"prime_values\": "<<np<<",\n  \"seven_prime_values\": "<<n7<<",\n  \"values\": [\n";
bool first=true;for(auto &[H,a]:rows){if(!first)cout<<",\n";first=false;cout<<"    {\"H\":"<<H<<",\"input_rows\":[";for(int i=0;i<(int)a.size();i++){if(i)cout<<",";cout<<a[i];}cout<<"],\"covers\":[";for(int i=0;i<3;i++){if(i)cout<<",";cout<<profiles[H][i];}cout<<"],\"least_divisor_or_one\":"<<factor(H)<<",\"seven_quotient_least_divisor_or_one\":"<<(H%7?0:factor(H/7))<<"}";}cout<<"\n  ]\n}\n";
}
