#include <bits/stdc++.h>
using namespace std;
struct Pair {int a,b;};
int main(int argc,char**argv){
  int n=argc>1?atoi(argv[1]):7;
  vector<Pair> ed;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)ed.push_back({i,j});
  int m=ed.size(),Q=1<<m,M=1<<n;
  vector<uint16_t> H(Q);
  int w[]={2,4,6,12};
  auto info=[&](int b,int *e){cout<<"n "<<n<<" base code "<<b<<" H "<<H[b]<<" edges";for(int k=0;k<4;k++)cout<<" ("<<ed[e[k]].a<<","<<ed[e[k]].b<<")";cout<<"\n";};
  for(int b=0;b<Q;b++){
    unsigned row[9]={};
    for(int i=0;i<m;i++){auto [a,z]=ed[i];if(b>>i&1)swap(a,z);row[a]|=1<<z;}
    unsigned short f[256][8]={};
    for(int i=0;i<n;i++)f[1<<i][i]=1;
    for(int mask=1;mask<M;mask++){
      if(!(mask&(mask-1)))continue;
      for(int u=mask;u;u&=u-1){
        auto i=__builtin_ctz(u); auto t=mask^(1<<i);int c=0;
        for(int v=t&row[i];v;v&=v-1)c+=f[t][__builtin_ctz(v)];
        f[mask][i]=c;
      }
    }
    int h=0; for(int i=0;i<n;i++)h+=f[M-1][i];H[b]=h;
  }
  if(n==7){
    for(int h:{65,67}){bool first=true;
    for(int b=0;b<Q && first;b++) if(H[b]==h)
      for(int i=0;i<m && first;i++)if(H[b^(1<<i)]==h+2)
        for(int j=0;j<m && first;j++)if(H[b^(1<<i)^(1<<j)]==h+4)
          for(int k=0;k<m && first;k++)if(H[b^(1<<i)^(1<<j)^(1<<k)]==h+6){
            cout<<"chain "<<h<<" "<<b<<" arcs ";
            for(auto e:{i,j,k})cout<<"("<<ed[e].a<<","<<ed[e].b<<")";
            cout<<"\n";first=false;
          }
    if(first)cout<<"no +2 chain from "<<h<<" through "<<h+6<<"\n";}
  }
  long long two=0,three=0,four=0;
  for(int b=0;b<Q;b++){
    vector<int> C[4];for(int i=0;i<m;i++){int dif=(int)H[b^(1<<i)]-H[b];for(int k=0;k<4;k++)if(dif==w[k])C[k].push_back(i);}
    for(int i:C[0])for(int j:C[1]){
      auto a=(1<<i)|(1<<j);
      if(H[b^a]!=H[b]+6)continue;
      two++;
      for(int k:C[2]){
        int c=1<<k;
        if(H[b^((1<<i)|c)]!=H[b]+8 || H[b^((1<<j)|c)]!=H[b]+10 || H[b^(a|c)]!=H[b]+12)continue;
        three++;
        for(int l:C[3]){
          int orig[]={i,j,k,l};bool good=true;
          for(int v=1;v<8;v++){
            int P=b^(1<<l),sum=12;
            for(int c=0;c<3;c++)if(v>>c&1){P^=1<<orig[c];sum+=w[c];}
            if(H[P]!=H[b]+sum){good=false;break;}
          }
          if(good){
            four++;if(four<=12)info(b,orig);
          }
        }
      }
    }
  }
  cout<<"checked "<<Q<<" labelled bases; two/three/four cubes "<<two<<" "<<three<<" "<<four<<"\n";
}
