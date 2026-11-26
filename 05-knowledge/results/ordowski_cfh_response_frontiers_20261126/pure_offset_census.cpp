
#include <algorithm>
#include <array>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <cstdint>
#include <numeric>
using namespace std;
using U=uint64_t;
const int MAX=9; int ix[MAX][MAX];
array<int,MAX> adj(int n,U w){
 array<int,MAX>a{};
 for(int j=1;j<n;j++)for(int i=0;i<j;i++)
  if(w>>ix[i][j]&1)a[j]|=1<<i;else a[i]|=1<<j;
 return a;
}
U canonical(int n,U w){
 auto a=adj(n,w);
 array<int,MAX> color{},p{};
 for(int i=0;i<n;i++)color[i]=__builtin_popcount(a[i]),p[i]=i;
 for(int step=0;step<n;step++){
  array<array<int,MAX+1>,MAX> s{};
  for(int i=0;i<n;i++){
   s[i][0]=color[i];for(int j=0;j<n;j++)if(a[i]>>j&1)s[i][1+color[j]]++;
  }
  sort(p.begin(),p.begin()+n,[&](int i,int j){return s[i]<s[j];});
  array<int,MAX> next{};
  for(int t=1;t<n;t++)next[p[t]]=next[p[t-1]]+(s[p[t]]!=s[p[t-1]]);
  const int old=*max_element(color.begin(),color.begin()+n), cl=next[p[n-1]];
  color=next;if(cl==n-1||(step>0&&old==cl))break;
 }
 vector<vector<int>> groups(n);for(int i=0;i<n;i++)groups[color[i]].push_back(i);
 groups.erase(remove_if(groups.begin(),groups.end(),[](const auto& v){return v.empty();}),groups.end());
 U best=~0ull;
 auto go=[&](auto&& f,int g,int j,U z)->void{
  if(g==(int)groups.size()){best=min(best,z);return;}
  auto& q=groups[g];
  do{
   int k=j;U mask=z;
   for(int u:q){p[k]=u;for(int i=0;i<k;i++)if(a[u]>>p[i]&1)mask|=1ull<<ix[i][k];k++;}
   if(mask<best)f(f,g+1,k,mask);
  }while(next_permutation(q.begin(),q.end()));
 };
 go(go,0,0,0);return best;
}
vector<U> extend(int n,const vector<U>& prev){
 unordered_set<U> C;
 for(U w:prev)for(U q=0;q<(1ull<<(n-1));q++)C.insert(canonical(n,w|(q<<ix[0][n-1])));
 vector<U> v(C.begin(),C.end());sort(v.begin(),v.end());return v;
}
int fwd[1<<MAX][MAX],revv[1<<MAX][MAX],P[1<<MAX];
struct Prof {int H;array<int,MAX>E,R,V,e;int Q;};
Prof profile(int n,U w,bool onlyH=false){
 auto A=adj(n,w); const int M=(1<<n)-1;
 for(int s=0;s<=M;s++)for(int j=0;j<n;j++)fwd[s][j]=0;
 for(int i=0;i<n;i++)fwd[1<<i][i]=1;
 P[0]=1;
 for(int s=1;s<=M;s++){
  int v=0;
  for(int i=0;i<n;i++)if(int z=fwd[s][i]){
   v+=z;for(int t=A[i]&~s;t;t&=t-1){int b=t&-t;fwd[s|b][__builtin_ctz(b)]+=z;}
  }
  P[s]=v;
 }
 Prof out{P[M],{},{},{},{},0};
 if(onlyH)return out;
 for(int s=0;s<=M;s++)out.Q+=P[s]*P[M^s];
 for(int s=0;s<=M;s++)for(int j=0;j<n;j++)revv[s][j]=0;
 for(int i=0;i<n;i++)revv[1<<i][i]=1;
 for(int s=1;s<=M;s++)for(int i=0;i<n;i++){
  if(int z=fwd[s][i])out.R[i]+=z*P[M^s];
  if(int z=revv[s][i]){
   out.V[i]+=z*P[M^s];
   for(int t=M&~(A[i]|s);t;t&=t-1){int b=t&-t;revv[s|b][__builtin_ctz(b)]+=z;}
  }
 }
 for(int j=0;j<n;j++){out.E[j]=revv[M][j];out.e[j]=fwd[M][j];}
 return out;
}
bool search(int n,const vector<U>& C){
 vector<pair<int,int>> arcs;
 constexpr int HH[7]={23,25,29,59,71,103,143};
 constexpr int RR[7]={24,28,31,87,72,129,173};
 for(int j=1;j<n;j++)for(int i=0;i<j;i++)arcs.push_back({i,j});
 uint64_t roots=0,marked=0,e1=0,opposite=0;
 int best=100000000, eligible=100000000;
 unordered_set<string> kinds;
 for(U w:C){
  Prof base=profile(n,w);
  bool yes=false;
  for(int r=0;r<n;r++)if(base.e[r]==1)yes=true,roots++;
  if(!yes)continue;
  for(int x=0;x<(int)arcs.size();x++){
   const U z=w^(1ull<<x);auto p=profile(n,z,true);if(p.H!=base.H)continue;
   p=profile(n,z);if(p.Q!=base.Q)continue;
   for(int r=0;r<n;r++)if(base.e[r]==1&&p.e[r]==1&&p.R[r]==base.R[r]){
    const int E=p.E[r]-base.E[r],V=p.V[r]-base.V[r];if(!E&&!V)continue;
    marked++;if(E*V<0)opposite++;
    const auto [u,v]=arcs[x];
    for(int k=0;k<7;k++){
     const int d=abs(HH[k]*E+RR[k]*V);
     if(d&&d<best){
      best=d;
      cout<<"witness "<<n<<" mask "<<w<<" root "<<r<<" arc "<<u<<","<<v<<" inner "<<base.H<<","<<base.R[r]<<","<<base.Q
          <<" dE,dV "<<E<<","<<V<<" outer "<<HH[k]<<","<<RR[k]<<" shift "<<d<<endl;
     }
     const long long A=HH[k]*base.H, B=HH[k]*base.R[r]+RR[k]*base.Q;
     if(d && d<eligible && gcd(A,B)==1){
      long long z=2*(B%A)%A;int bit=1;
      while(z!=B%A&&z!=2&&z!=A-2){z=2*z%A;bit++;}
      if(z==2||z==A-2){
       eligible=d;
       cout<<"unit-ready "<<n<<" mask "<<w<<" root "<<r<<" arc "<<u<<","<<v<<" inner "<<base.H<<","<<base.R[r]<<","<<base.Q
         <<" dE,dV "<<E<<","<<V<<" outer "<<HH[k]<<","<<RR[k]<<" offset "<<d<<" spoke "<<bit<<","<<(z==2?2:-2)<<endl;
      }
     }
     if(d==2)e1++;
    }
   }
  }
 }
 cout<<"n "<<n<<" classes "<<C.size()<<" unique ending roots "<<roots<<" pure offset arcs "<<marked<<" opposing terms "<<opposite<<" best "<<best<<" unit-ready best "<<eligible<<" unit "<<e1<<endl;
 return false;
}
int main(int argc,char**argv){
 int n=argc>1?atoi(argv[1]):8;
 if(n<2||n>MAX)return 1;
 for(int j=1;j<MAX;j++)for(int i=0;i<j;i++)ix[i][j]=j*(j-1)/2+i;
 vector<U> C{0};for(int t=2;t<=n;t++){C=extend(t,C);if(search(t,C))break;}
}
