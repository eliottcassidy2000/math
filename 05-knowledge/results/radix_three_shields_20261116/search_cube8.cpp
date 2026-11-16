#include <bits/stdc++.h>
using namespace std;
int main(int argc,char**argv){
  if(argc!=2) return 2;
  vector<pair<int,int>> edges;for(int i=0;i<8;i++)for(int j=i+1;j<8;j++)edges.push_back({i,j});
  auto count=[&](unsigned bits){
    unsigned rows[8]={};
    for(int j=0;j<28;j++){auto [a,b]=edges[j];if(bits>>j&1)swap(a,b);rows[a]|=1u<<b;}
    int p[256][8]={};
    for(int i=0;i<8;i++)p[1u<<i][i]=1;
    for(int mask=1;mask<256;mask++){
      if(!(mask&(mask-1)))continue;
      for(int t=mask;t;t&=t-1){
        int i=__builtin_ctz(t),old=mask^(1<<i);
        for(int s=old&rows[i];s;s&=s-1)p[mask][i]+=p[old][__builtin_ctz(s)];
      }
    }
    return accumulate(p[255],p[255]+8,0);
  };
  ifstream in(argv[1]);string line;getline(in,line);
  int w[]={2,4,6,12};
  long long n=0, two=0, three=0, four=0;
  bool found[2]={};
  while(getline(in,line)){
    if(line.empty())continue;
    unsigned long long rows=stoull(line);
    unsigned b=0;
    for(int i=0;i<28;i++){auto [x,y]=edges[i];if(!(rows>>(8*x+y)&1))b|=1u<<i;}
    int H=count(b);n++;
    if((H==65||H==67)&&!found[(H-65)/2]){
      for(unsigned i=0;i<28&&!found[(H-65)/2];i++)
        if(count(b^(1u<<i))==H+2)
          for(unsigned j=0;j<28&&!found[(H-65)/2];j++)
            if(count(b^(1u<<i)^(1u<<j))==H+4)
              for(unsigned k=0;k<28&&!found[(H-65)/2];k++)
                if(count(b^(1u<<i)^(1u<<j)^(1u<<k))==H+6){
                  found[(H-65)/2]=true;
                  cout<<"chain "<<H<<" "<<b<<" arcs ";
                  for(auto z:{i,j,k})cout<<"("<<edges[z].first<<","<<edges[z].second<<")";
                  cout<<"\n";
                }
    }
    vector<unsigned>C[4];for(unsigned i=0;i<28;i++){
      int d=count(b^(1u<<i))-H;
      for(int k=0;k<4;k++)if(d==w[k])C[k].push_back(1u<<i);
    }
    for(auto i:C[0])for(auto j:C[1]){
      if(count(b^i^j)!=H+6)continue;two++;
      for(auto k:C[2]){
        if(count(b^i^k)!=H+8||count(b^j^k)!=H+10||count(b^i^j^k)!=H+12)continue;
        three++; if(three<=4)cout<<"three "<<b<<" H "<<H<<" masks "<<i<<" "<<j<<" "<<k<<"\n";
        for(auto l:C[3]){
          bool good=true;unsigned E[]={i,j,k};
          for(int z=0;z<8;z++){
            unsigned u=b^l;int total=H+12;
            for(int f=0;f<3;f++)if(z>>f&1){u^=E[f];total+=w[f];}
            if(count(u)!=total){good=false;break;}
          }
          if(good){
            four++; if(four<=4)cout<<"four "<<b<<" H "<<H<<" masks "<<i<<" "<<j<<" "<<k<<" "<<l<<"\n";
          }
        }
      }
    }
  }
  for(int i=0;i<2;i++)if(!found[i])cout<<"no +2 chain from "<<65+2*i<<" through "<<71+2*i<<"\n";
  cout<<"n8 representatives "<<n<<"; two/three/four cubes "<<two<<" "<<three<<" "<<four<<"\n";
}
