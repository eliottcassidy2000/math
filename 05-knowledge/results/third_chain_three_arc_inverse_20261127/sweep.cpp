#include <bits/stdc++.h>
using namespace std;using I=__int128_t;
struct T{int k;int e[5];};
unsigned __int128 v[5][100]; long long tar[39];int modA=65521,modB=65519;
int p[2][5][100],gap[5],maxn,ct=0;
long long visits=0,matches[60]={};unordered_map<int,vector<int>> candidates;
vector<T> terms; vector<pair<int,int>> B; int m, from[6];
void rec(int g){
 if(g==m-1){
 ++visits;
 int r=0;
 for(auto &t:terms){int u=t.k;for(int j=0;j<m-1;j++)u=int((long long)u*p[0][t.e[j]][gap[j]]%modA);r+=u;}r=(r%modA+modA)%modA;
 auto a=candidates.find(r); if(a==candidates.end())return;
 int s=0;for(auto &t:terms){int u=t.k;for(int j=0;j<m-1;j++)u=int((long long)u*p[1][t.e[j]][gap[j]]%modB);s+=u;}s=(s%modB+modB)%modB;
 for(int c:a->second)if(s==tar[c]%modB){
  unsigned __int128 n=0;for(auto &t:terms){unsigned __int128 u=t.k;for(int j=0;j<m-1;j++)u*=v[t.e[j]][gap[j]];n+=u;}
  if(n==tar[c]){if(matches[c]++==0){cout<<c<<" witness ";for(auto e:B)cout<<from[e.first]<<","<<from[e.second]<<" ";cout<<" order "<<from[m-1]+1<<endl;}}
 }return;
 }
 int l=0,h=maxn;bool active=false;
 for(auto [i,j]:B)if(i<=g&&g<j){active=true;h=min(h,maxn-(from[g]-from[i])+(g-i));if(i==g&&j==g+1)l=1;}
 if(!active)h=0;
 for(int x=l;x<=h;x++){gap[g]=x;from[g+1]=from[g]+x+1;rec(g+1);}
}
int main(int argc,char**argv){
 if(argc!=2)return 2;int stop=atoi(argv[1]);if(stop<0||stop>23)return 2;tar[0]=7;
 for(int t=0;t<=stop;t++){if(t)tar[t]=tar[t-1]*3;candidates[tar[t]%modA].push_back(t);}
 I u=tar[stop]-1; while(u>1)u>>=1,maxn++;
 // The cap is at most 39. Each of at most 720 exceptional orders has
 // at most 2^(3*cap) slot assignments, so true counts fit unsigned 128 bits.
 // Unsigned arithmetic on signed monomials computes the exact result mod 2^128.
 for(int j=0;j<=4;j++){v[j][0]=p[0][j][0]=p[1][j][0]=1;for(int h=1;h<=maxn;h++){v[j][h]=v[j][h-1]*j;p[0][j][h]=p[0][j][h-1]*j%modA;p[1][j][h]=p[1][j][h-1]*j%modB;}}
 int nb,k;cin>>nb;while(nb--){cin>>m>>k;B.resize(3);for(auto &e:B)cin>>e.first>>e.second;terms.resize(k);for(auto &t:terms){cin>>t.k;for(int g=0;g<m-1;g++)cin>>t.e[g];}rec(0);}
 cout<<"visits "<<visits<<" common cap "<<maxn<<endl;
 for(int t=0;t<=stop;t++)cout<<t<<":"<<matches[t]<<" ";cout<<endl;
}
