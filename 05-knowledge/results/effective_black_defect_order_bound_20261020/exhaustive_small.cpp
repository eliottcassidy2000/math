#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>
using namespace std;

int main() {
  cout << "Independent labelled-tournament scan; direct permutations and reachability\n";
  for(int n=3;n<=6;n++) {
    int e=n*(n-1)/2, full=(1<<n)-1;
    vector<vector<int>> ps;
    vector<int> p(n); iota(p.begin(),p.end(),0);
    do {ps.push_back(p);} while(next_permutation(p.begin(),p.end()));
    uint64_t strong_count=0, cycles=0, q1cycles=0,q2pluscycles=0;
    map<int,uint64_t> defects, qhist;
    for(int mask=0;mask<(1<<e);mask++) {
      bool a[6][6]={}; int out[6]={};
      for(int i=0,k=0;i<n;i++) for(int j=i+1;j<n;j++,k++) {
        a[i][j]=(mask>>k)&1; a[j][i]=!a[i][j];
        if(a[i][j])out[i]|=1<<j; else out[j]|=1<<i;
      }
      auto strong=[&](int set) {
        if(!set)return false;
        int first=__builtin_ctz((unsigned)set);
        for(int rev=0;rev<2;rev++) {
          int seen=1<<first;
          while(true) {
            int next=seen;
            for(int u=0;u<n;u++)if(seen>>u&1)
              next|=(rev ? (full^(out[u]|(1<<u))) : out[u])&set;
            if(next==seen)break; seen=next;
          }
          if(seen!=set)return false;
        }
        return true;
      };
      if(!strong(full))continue;
      strong_count++;
      int kappa=n;
      for(int keep=1;keep<=full;keep++)if(!strong(keep))
        kappa=min(kappa,n-__builtin_popcount((unsigned)keep));
      int aut=0,H=0,bluewords=0;
      vector<vector<int>> hcycles;
      for(const auto &p:ps) {
        bool autom=true;
        for(int i=0;i<n && autom;i++)for(int j=i+1;j<n;j++)
          if(a[i][j]!=a[p[i]][p[j]]){autom=false;break;}
        aut+=autom;
        bool path=true;
        for(int i=0;i<n-1;i++)if(!a[p[i]][p[i+1]]){path=false;break;}
        if(!path)continue;
        H++;
        bool blue=true;
        for(int i=0;i<n && blue;i++)for(int j=i+1;j<n;j++)
          if(a[p[i]][p[j]]==a[p[n-1-i]][p[n-1-j]]){blue=false;break;}
        bluewords+=blue;
        if(p[0]==0 && a[p[n-1]][0])hcycles.push_back(p);
      }
      assert(aut%2==1 && H%aut==0 && bluewords%aut==0);
      int Delta=(H-bluewords)/aut;
      defects[Delta]++;
      int s=Delta+2,r=1; while(r<4*s*(Delta+1))r*=2;
      assert(n<=2*s*r && 2*s*r<16*s*s*(Delta+1));
      for(const auto &p:hcycles) {
        cycles++;
        int m=0;
        for(int shift=0;shift<n;shift++) {
          bool ok=true;
          for(int i=0;i<n && ok;i++)for(int j=i+1;j<n;j++)
            if(a[p[i]][p[j]]!=a[p[(i+shift)%n]][p[(j+shift)%n]]){ok=false;break;}
          m+=ok;
        }
        assert(m>0 && m%2==1 && n%m==0);
        int q=n/m; qhist[q]++;
        int bluecuts=0;
        for(int shift=0;shift<q;shift++) {
          bool ok=true;
          for(int i=0;i<n && ok;i++)for(int j=i+1;j<n;j++)
            if(a[p[(shift+i)%n]][p[(shift+j)%n]]==a[p[(shift+n-1-i)%n]][p[(shift+n-1-j)%n]]) {ok=false;break;}
          bluecuts+=ok;
        }
        assert(bluecuts==0 || bluecuts==gcd(2,q));
        assert(q-bluecuts<=Delta && q<=Delta+2);
        for(int i=0;i<n;i++)assert(__builtin_popcount((unsigned)out[i])>=(m-1)/2);
        if(q==1){q1cycles++; assert(kappa>=(m+2)/3);}
        else {q2pluscycles++; assert(kappa>=(m+1)/2);}
      }
    }
    cout<<"n="<<n<<" strong="<<strong_count<<" cycles="<<cycles<<" q1="<<q1cycles<<" q_ge_2="<<q2pluscycles<<"\n";
    cout<<"defect_histogram:"; for(auto [d,c]:defects)cout<<" "<<d<<":"<<c; cout<<"\n";
    cout<<"cycle_q_histogram:"; for(auto [q,c]:qhist)cout<<" "<<q<<":"<<c;cout<<"\n";
  }
  cout<<"PASS: all assertions, n=3,...,6\n";
}
