// A focused birth-splice route: restore a source, then enlarge the same TT module.
#define main previous_search_main
#include "../adaptive_four_returns_20261103/lift_search.cpp"
#undef main
int main(int argc,char**argv){
 if(argc<3)return 2;ifstream in(argv[1]);int n,port;in>>n>>port;T t(n);for(auto&v:t)in>>v;
 int limit=stoi(argv[2]);I h=dp(t).h.back(),h2=dp(clone(t,port)).h.back(),B=2*h-h2;
 I base=(9*h-3*h2)/2;
 if((9*h-3*h2)%2){cerr<<"invalid source parity\n";return 4;}
 if(base%2==0){cerr<<"excluded by intermediate Hamilton parity "<<base<<'\n';return 0;}
 struct C{I h;T p,t;};vector<C> first;
 for(auto&p:paths(t)){T u=native(t,p);first.push_back({dp(u).h.back(),p,u});}
 sort(first.begin(),first.end(),[](auto&a,auto&b){return a.h!=b.h?a.h<b.h:a.p<b.p;});
 auto frame=[&](const T&p){int i=find(p.begin(),p.end(),port)-p.begin();return pair<int,int>{i?p[i-1]:-1,i+1<n?p[i+1]:-1};};
 bool found[2]={false,false};I examined=0,scalar=0;
 for(int rank=0;rank<min(limit,(int)first.size());rank++){
  auto&a=first[rank];cerr<<"rank "<<rank<<'\n';
  for(auto&q:paths(a.t))if(frame(a.p)==frame(q)){
   T u=native(a.t,q);I v=dp(u).h.back();examined++;
   int lane=v==base+2?0:v==base+4?1:-1;if(lane<0||found[lane]||!strong(u))continue;scalar++;
   I z=dp(clone(u,port)).h.back();if(z!=3*h+2+2*lane)continue;
   // Append the newborn after the substituted block; all exterior relations are uniform.
   int mask=1<<port;for(int i=0;i<n;i++)if(u[i]>>port&1)mask|=1<<i;
   T w=birth(u,mask);cout<<"{\"source_H\":"<<h<<",\"port\":"<<port<<",\"offset\":"<<2+2*lane<<",\"B\":"<<B<<",\"intermediate_H\":"<<v<<",\"rank\":"<<rank<<",\"examined\":"<<examined<<",\"first_path\":";print(a.p);cout<<",\"second_path\":";print(q);cout<<",\"mask\":"<<mask<<",\"target\":";print(w);cout<<"}\n";cout.flush();found[lane]=true;
   if(found[0]&&found[1])return 0;
  }
 }
 cerr<<"bounded end examined "<<examined<<" intermediate hits "<<scalar<<" found "<<found[0]<<' '<<found[1]<<'\n';return 0;
}
