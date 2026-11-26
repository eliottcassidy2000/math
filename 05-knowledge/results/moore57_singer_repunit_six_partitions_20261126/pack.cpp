#include <array>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdint>
using namespace std;
using U=uint64_t; using R=array<int,6>;
struct Pack {
  static constexpr U ALL=(1ull<<54)-1;
  U data[6][6][54]{}, live[6][6][54]{};
  mt19937_64 rng;
  vector<R> rs;
  int visit=0,limit=3000;
  Pack(string file,int seed):rng(seed){
    ifstream f(file);
    for(int i=0;i<6;i++)for(int j=i+1;j<6;j++)for(int x=0;x<54;x++){
      string word;f>>word; U y=stoull(word,0,16);
      data[i][j][x]=y;while(y){int a=__builtin_ctzll(y);data[j][i][a]|=1ull<<x;y&=y-1;}
    }
    for(int i=0;i<6;i++)for(int j=0;j<6;j++)for(int x=0;x<54;x++)live[i][j][x]=data[i][j][x];
  }
  int rnd(int n){return rng()%n;}
  int randbit(U a){int n=rnd(__builtin_popcountll(a));while(n--)a&=a-1;return __builtin_ctzll(a);}
  void toggle(R a,int d){
    for(int i=0;i<6;i++)for(int j=i+1;j<6;j++){
      U v=1ull<<a[j],w=1ull<<a[i];
      if(d==1){live[i][j][a[i]]&=~v;live[j][i][a[j]]&=~w;}
      else{live[i][j][a[i]]|=v;live[j][i][a[j]]|=w;}
    }
  }
  bool dfs(R& a,array<U,6> d,int left){
    if(!left)return true;
    if(++visit>limit)return false;
    int k=-1,minn=55;
    for(int j=0;j<6;j++)if((left>>j)&1){
      int p=__builtin_popcountll(d[j]);if(p<minn){minn=p;k=j;}
    }
    if(!minn)return false;
    left^=1<<k;
    U u=d[k];
    while(u){
      int x=randbit(u);u^=1ull<<x;a[k]=x;
      auto dd=d;
      for(int j=0;j<6;j++)if((left>>j)&1)dd[j]&=live[k][j][x];
      if(dfs(a,dd,left))return true;
      if(visit>limit)break;
    }
    return false;
  }
  bool start(int i,int x,R& a){
    int left=63^(1<<i);a[i]=x;
    array<U,6> d;
    for(int j=0;j<6;j++)d[j]=live[i][j][x];
    visit=0;return dfs(a,d,left);
  }
  int sweep(){
    vector<int> labels(324);for(int k=0;k<324;k++)labels[k]=k;
    shuffle(labels.begin(),labels.end(),rng);
    stable_sort(labels.begin(),labels.end(),[&](int a,int b){
      int i=a/54,x=a%54,j=b/54,y=b%54;
      return __builtin_popcountll(live[i][(i+1)%6][x])<__builtin_popcountll(live[j][(j+1)%6][y]);
    });
    for(int k:labels){
      R a; if(start(k/54,k%54,a)){rs.push_back(a);toggle(a,1);return 1;}
    }
    return 0;
  }
  void report(const string& out){
    ofstream f(out);
    for(const auto& a:rs){for(int x:a)f<<x<<" ";f<<"\n";}
  }
  void migrate(){
    for(int it=0;it<256;it++){
      int j=rnd(rs.size()),i=rnd(6); R a=rs[j]; U able=ALL;
      for(int k=0;k<6;k++)if(k!=i)able &= live[k][i][a[k]];
      if(able){
        toggle(a,-1);a[i]=randbit(able);
        rs[j]=a;toggle(a,1);
      }
    }
  }
  void load(const string &out){
    ifstream f(out); R a;
    while(f>>a[0]){
      for(int j=1;j<6;j++)f>>a[j];
      for(int i=0;i<6;i++)for(int j=i+1;j<6;j++)
        if(!((live[i][j][a[i]]>>a[j])&1))throw runtime_error("bad seed pack");
      rs.push_back(a);toggle(a,1);
    }
  }
  void greedy(){
    while(sweep()){}
  }
  void go(const string& output,int iterations){
    greedy();int best=rs.size();report(output);
    cerr<<"initial "<<best<<"\n";
    for(int k=0;k<iterations;k++){
      migrate();
      int orig=rs.size(), del=1+rnd(k%9==8?24:7);
      vector<R> trash;
      for(int t=0;t<del;t++){
        int j=rnd(rs.size());auto a=rs[j];trash.push_back(a);toggle(a,-1);
        swap(rs[j],rs.back());rs.pop_back();
      }
      int mid=rs.size();
      greedy();
      if(rs.size()<orig && (int(rs.size())+2<orig||rnd(128)!=0)){
        while((int)rs.size()>mid){toggle(rs.back(),-1);rs.pop_back();}
        for(const auto& a:trash){rs.push_back(a);toggle(a,1);}
      }
      if((int)rs.size()>best){best=rs.size();cerr<<k<<" "<<best<<"\n";report(output);}
    }
  }
};
int main(int argc,char**argv){
  if(argc<5||argc>6)return 1;
  Pack p(argv[1],stoi(argv[3]));if(argc==6)p.load(argv[5]);p.go(argv[2],stoi(argv[4]));
}
