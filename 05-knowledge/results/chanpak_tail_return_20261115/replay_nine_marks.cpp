#define main ordinary_nine_main
#include "replay_nine_census.cpp"
#undef main
#include <tuple>
#include <set>
using Key=tuple<int,int,int>;
int main(){
  unordered_set<U> old={0};
  for(n=2;n<=8;n++){
    unordered_set<U> next;
    for(U c:old){
      Row r=decode(c,n-1);
      for(int b=0;b<(1<<(n-1));b++){
        Row s=r;s[n-1]=b;
        for(int i=0;i<n-1;i++)if(!((b>>i)&1))s[i]|=1<<(n-1);
        next.insert(canonical(s,{(1<<n)-1}));
      }
    }
    old=move(next);
  }
  if(old.size()!=6880) return 1;
  const int F=255;
  map<Key,set<int>> idx;
  map<Key,map<int,pair<Row,int>>> wit;
  map<int,int> seen;
  map<pair<int,int>,map<int,pair<Row,int>>> normalized;
  struct Record{int h,S,d,H;U core;int mark;bool eligible;};
  vector<Record> ALL;
  for(U c:old){
    Row r=decode(c,8);
    int paths[256][8]={}, h[256]={1};
    for(int m=1;m<=F;m++) for(int b=m;b;b&=b-1) {
      int v=__builtin_ctz(unsigned(b));int s=m^(1<<v);
      if(!s)paths[m][v]=1;
      else for(int t=s&r[v];t;t&=t-1) paths[m][v]+=paths[s][__builtin_ctz(unsigned(t))];
      h[m]+=paths[m][v];
    }
    int SS[8]={},dd[8]={},P[8][8]={};
    for(int a=0;a<8;a++){
      int way[256][8]={};
      way[1<<a][a]=1;
      for(int m=1;m<=F;m++)if((m>>a)&1){
        int l=__builtin_popcount(unsigned(m));
        int p=F^m;
        SS[a]+=paths[m][a]*h[p];
        dd[a]+=(l&1 ? 1:-1)*paths[m][a]*h[p];
        for(int b=m;b;b&=b-1){
          int v=__builtin_ctz(unsigned(b));
          int x=way[m][v];
          if(!x)continue;
          if(!(l&1))P[a][v]+=h[p]*x;
          for(int t=p&r[v];t;t&=t-1){
            int j=__builtin_ctz(unsigned(t));way[m|(1<<j)][j]+=x;
          }
        }
      }
    }
    for(int b=0;b<=F;b++){
      int H=h[F],S=0,d=0;
      for(int a=b;a;a&=a-1){
        int i=__builtin_ctz(unsigned(a));
        S+=SS[i];d+=dd[i];
        for(int t=F^b;t;t&=t-1)H+=2*P[i][__builtin_ctz(unsigned(t))];
      }
      if(c%100==0 && b%16==0) {
        Row g=r;g[8]=b;
        for(int i=0;i<8;i++)if(!((b>>i)&1))g[i]|=1<<8;
        if(H!=ham(g))return 2;
      }
      Key z={h[F],S,d};
      int reach=b;
      for(int i=0;i<8;i++)for(int a=reach;a;a&=a-1)reach|=r[__builtin_ctz(unsigned(a))];
      ALL.push_back({h[F],S,d,H,c,b,reach==F});
      int a=S+h[F], lam=a/(a&-a);
      if(H>=1255&&H<=1271){
        pair<Row,int> v={r,b};
        auto &M=normalized[{h[F],lam}];
        if(!M.count(H)||v<M[H])M[H]=v;
      }
      if(H>=1255&&H<=1271) {
        idx[z].insert(H);seen[H]++;
      }
      if(H>=1255&&H<=1271) {
        pair<Row,int> v={r,b};
        if(!wit[z].count(H)||v<wit[z][H])wit[z][H]=v;
      }
    }
  }
  cout<<"Surgery kernel checked against direct path DP on proportional controls"<<endl;
  cout<<"Target marked extensions (nonquotiented):";
  for(auto [H,k]:seen)cout<<" "<<H<<":"<<k;
  cout<<endl;
  int best=0,contig=0,nn=0;
  Key kbest={};
  for(auto &[z,si]:idx){
    best=max(best,int(si.size()));
    for(int H:si)if(!si.count(H-2)){
      int len=1;while(si.count(H+2*len))len++;
      if(len>contig){contig=len;kbest=z;nn=H;}
    }
  }
  cout<<"Common right-state profiles "<<idx.size()<<" maximum distinct lane values "<<best<<" longest adjacent run "<<contig<<endl;
  auto [h,S,d]=kbest;
  cout<<"Run state "<<h<<" "<<S<<" "<<d<<":";
  for(int i=0;i<contig;i++)cout<<" "<<nn+2*i;
  cout<<endl;
  for(int i=0;i<contig;i++){
    auto [r,b]=wit[kbest][nn+2*i];
    cout<<"  "<<nn+2*i<<" rows";
    for(int j=0;j<8;j++)cout<<" "<<r[j];
    cout<<" mark "<<b<<endl;
  }
  int nr=0,na=0;
  pair<int,int> nk={};
  for(auto &[z,si]:normalized)for(auto &[H,w]:si)if(!si.count(H-2)){
    int len=1;while(si.count(H+2*len))len++;
    if(len>nr){nr=len;na=H;nk=z;}
  }
  cout<<"Longest first-rung run after neutral-one alignment "<<nr<<" normalized h,oddpart(S+h) "<<nk.first<<" "<<nk.second<<" contiguous from "<<na<<endl;
  for(int i=0;i<nr;i++){
    auto [r,b]=normalized[nk][na+2*i];
    cout<<"  "<<na+2*i<<" rows";
    for(int j=0;j<8;j++)cout<<" "<<r[j];
    cout<<" mark "<<b<<endl;
  }
  sort(ALL.begin(),ALL.end(),[](const Record&a,const Record&b){
    return tie(a.h,a.S,a.d,a.H,a.core,a.mark)<tie(b.h,b.S,b.d,b.H,b.core,b.mark);
  });
  int br=0,curr=0,prev=-2;
  Key key={};
  Record brstart={},runstart={};
  for(const auto &a:ALL){
    Key z={a.h,a.S,a.d};
    if(z==key&&a.H==prev)continue;
    if(z==key&&a.H==prev+2)curr++;
    else {curr=1;runstart=a;}
    prev=a.H;key=z;
    if(curr>br){br=curr;brstart=runstart;}
  }
  cout<<"Maximum exact global right-state run at order 9 "<<br<<" state "<<brstart.h<<" "<<brstart.S<<" "<<brstart.d<<" contiguous from "<<brstart.H<<endl;
  for(int i=0;i<br;i++){
    auto it=find_if(ALL.begin(),ALL.end(),[&](const Record &a){
      return a.h==brstart.h&&a.S==brstart.S&&a.d==brstart.d&&a.H==brstart.H+2*i;
    });
    Row r=decode(it->core,8);
    cout<<"  "<<it->H<<" rows";
    for(int j=0;j<8;j++)cout<<" "<<r[j];
    cout<<" mark "<<it->mark<<endl;
  }
  for(auto &a:ALL){int x=a.S+a.h;a.d=x/(x&-x);}
  sort(ALL.begin(),ALL.end(),[](const Record&a,const Record&b){
    return tie(a.h,a.d,a.H,a.S,a.core,a.mark)<tie(b.h,b.d,b.H,b.S,b.core,b.mark);
  });
  br=0;curr=0;prev=-2;key={};
  for(const auto &a:ALL){
    Key z={a.h,0,a.d};
    if(z==key&&a.H==prev)continue;
    if(z==key&&a.H==prev+2)curr++;
    else {curr=1;runstart=a;}
    prev=a.H;key=z;
    if(curr>br){br=curr;brstart=runstart;}
  }
  cout<<"Maximum global neutral-one alignable run at order 9 "<<br<<" h,oddpart(S+h) "<<brstart.h<<" "<<brstart.d<<" contiguous from "<<brstart.H<<endl;
  for(int i=0;i<br;i++){
    auto it=find_if(ALL.begin(),ALL.end(),[&](const Record &a){
      return a.h==brstart.h&&a.d==brstart.d&&a.H==brstart.H+2*i;
    });
    Row r=decode(it->core,8);
    cout<<"  "<<it->H<<" rows";
    for(int j=0;j<8;j++)cout<<" "<<r[j];
    cout<<" mark "<<it->mark<<" S+h "<<it->S+it->h<<endl;
  }

  br=0;curr=0;prev=-2;key={};
  for(const auto &a:ALL){
    if(!a.eligible)continue;
    Key z={a.h,0,a.d};
    if(z==key&&a.H==prev)continue;
    if(z==key&&a.H==prev+2)curr++;
    else {curr=1;runstart=a;}
    prev=a.H;key=z;
    if(curr>br){br=curr;brstart=runstart;}
  }
  cout<<"Maximum strong-completable neutral-one alignable run at order 9 "<<br<<" h,oddpart(S+h) "<<brstart.h<<" "<<brstart.d<<" contiguous from "<<brstart.H<<endl;
  for(int i=0;i<br;i++){
    auto it=find_if(ALL.begin(),ALL.end(),[&](const Record &a){
      return a.eligible&&a.h==brstart.h&&a.d==brstart.d&&a.H==brstart.H+2*i;
    });
    Row r=decode(it->core,8);
    cout<<"  "<<it->H<<" rows";
    for(int j=0;j<8;j++)cout<<" "<<r[j];
    cout<<" mark "<<it->mark<<" S+h "<<it->S+it->h<<endl;
  }
  return 0;

}
