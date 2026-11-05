#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <fstream>
#include <vector>
using namespace std;
int main(int argc, char** argv){
  map<pair<int,int>,long long> equal,boundary;
  set<array<int,8>> signatures;
  long long graphs=0,roots=0,zeros=0,ineqs=0;
  for(int n=1;n<=6;n++){
    int M=n*(n-1)/2; long long ng=0,nz=0;
    for(int code=0;code<(1<<M);code++){
      vector<int> adj(n),e(n),R(n),d(n),L(n),B(n),start(n); int h=0,k=0;
      for(int a=0;a<n;a++)for(int b=a+1;b<n;b++,k++)
        if(code&(1<<k))adj[a]|=1<<b;else adj[b]|=1<<a;
      vector<int> p(n);iota(p.begin(),p.end(),0);
      do{
        int bad=0,at=-1;
        for(int j=0;j<n-1;j++)if(!(adj[p[j]]&(1<<p[j+1]))){bad++;at=j;}
        if(bad==0){h++;e[p.back()]++;start[p[0]]++;for(int r:p){R[r]++;L[r]++;}}
        else if(bad==1){R[p[at+1]]++;L[p[at]]++;}
        int pb=bad - (n>1 && !(adj[p[n-2]]&(1<<p.back())));
        if(pb==0){d[p.back()]++;B[p.back()]+=n;}
        else if(pb==1)B[p.back()]++;
      }while(next_permutation(p.begin(),p.end()));
      for(int r=0;r<n;r++){
        roots++;
        signatures.insert({h,start[r],L[r],e[r],R[r],B[r],n==1,d[r]});
        if(!(L[r]>=h && h>=start[r] && R[r]>=h && B[r]>=1)){cerr<<"gadget bound fail\n";return 5;}
        if(h==1){
          int j=n-1-__builtin_popcount((unsigned)adj[r]);
          if(start[r]!=(j==0) || L[r]!=(1<<j) || e[r]!=(j==n-1) || R[r]!=(1<<(n-1-j)) || B[r]!=(1<<(n-1))){cerr<<"transitive formula fail\n";return 6;}
        }
        bool slack=(R[r]+e[r]-h-d[r]==0), small=(__builtin_popcount((unsigned)adj[r])<=1);
        if(slack!=small){cerr<<"zero slack mismatch "<<n<<' '<<code<<' '<<r<<'\n';return 1;}
        if(!(0<=e[r] && e[r]<=d[r] && d[r]<=h && h<=R[r] && R[r]+e[r]>=h+d[r])){cerr<<"inequality fail\n";return 2;}
        ineqs++;zeros+=slack;nz+=slack;
        if(R[r]==h){
          if(e[r]!=d[r]){cerr<<"equal response endpoint fail\n";return 3;}
          if(h>d[r]){
            if(h<3*d[r]){cerr<<"ratio fail "<<n<<' '<<code<<' '<<r<<'\n';return 4;}
            equal[{h,d[r]}]++;
          }else boundary[{h,d[r]}]++;
        }
      }
      graphs++;ng++;
    }
    cerr<<"n="<<n<<" graphs="<<ng<<" zero_slack_roots="<<nz<<"\n";
  }
  ofstream sig(argc>1 ? argv[1] : "bank6.txt");
  if(!sig){cerr<<"cannot write signature file\n";return 7;}
  for(auto a:signatures){for(int x:a)sig<<x<<' ';sig<<'\n';}
  cout<<"signatures "<<signatures.size()<<"\n";
  cout<<"graphs "<<graphs<<" roots "<<roots<<" zero_slack_roots "<<zeros<<" all_physical_inequalities "<<ineqs<<"\n";
  cout<<"nontrivial_equal_response\n";
  for(auto [p,m]:equal)cout<<p.first<<' '<<p.second<<' '<<m<<'\n';
  cout<<"boundary_equal_response\n";
  for(auto [p,m]:boundary)cout<<p.first<<' '<<p.second<<' '<<m<<'\n';
}
