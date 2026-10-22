#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>
using namespace std;
int main(){
 for(int n=1;n<=5;n++){
  int E=n*(n-1)/2, ids[5][5]={},lo[10],hi[10],k=0;
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){ids[i][j]=ids[j][i]=k;lo[k]=i;hi[k++]=j;}
  vector<vector<int>> ps; vector<int> p(n);iota(p.begin(),p.end(),0);do{ps.push_back(p);}while(next_permutation(p.begin(),p.end()));
  map<int,long long> hist, relaxedhist, maxperQhist; int overall=0,example=0;
  for(int q=0;q<(1<<E);q++){
   int ctrl[10],relaxed=0,qmax=0;
   for(int e=0;e<E;e++)ctrl[e]=(q>>e&1)?hi[e]:lo[e];
   for(auto &p:ps){int pinned=1<<p[n-1];for(int i=0;i+1<n;i++)pinned|=1<<ctrl[ids[p[i]][p[i+1]]];relaxed+=1<<(n-__builtin_popcount((unsigned)pinned));}
   relaxedhist[relaxed]++;
   for(int lab=0;lab<(1<<E);lab++){
    int F=0;
    for(auto &p:ps){int pinned=1<<p[n-1],plus=pinned;bool ok=true;
     for(int i=0;i+1<n;i++){int a=p[i],b=p[i+1],e=ids[a][b],v=ctrl[e],bit=1<<v;bool sign=((lab>>e)&1)^(a>b);if((pinned&bit)&&(((plus&bit)!=0)!=sign)){ok=false;break;}pinned|=bit;if(sign)plus|=bit;}
     if(ok)F+=1<<(n-__builtin_popcount((unsigned)pinned));
    }
    hist[F]++;qmax=max(qmax,F);if(F>overall){overall=F;example=q|(lab<<E);}
   }
   maxperQhist[qmax]++;
  }
  cout<<"n="<<n<<" max="<<overall<<" example="<<example<<"\n actual histogram:";for(auto [f,c]:hist)cout<<" "<<f<<":"<<c;cout<<"\n relaxed histogram:";for(auto[f,c]:relaxedhist)cout<<" "<<f<<":"<<c;cout<<"\n max-per-Q histogram:";for(auto[f,c]:maxperQhist)cout<<" "<<f<<":"<<c;cout<<"\n";
 }
}
