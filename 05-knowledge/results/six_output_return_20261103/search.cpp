// Positive searches at the six output heads, preserving the original clone port.
#define main previous_search_main
#include "../adaptive_four_returns_20261103/lift_search.cpp"
#undef main
int main(int argc,char**argv){
 if(argc<3)return 2;ifstream in(argv[1]);int n,port;in>>n>>port;T t(n);for(auto&v:t)in>>v;
 int limit=stoi(argv[2]),start=argc>3?stoi(argv[3]):0;bool widened=argc>4;I h=dp(t).h.back(),h2=dp(clone(t,port)).h.back();
 bool incoming=popcount((unsigned)t[port])==n-2;
 if(!incoming&&popcount((unsigned)t[port])!=1)return 3;
 struct C{I h;T p,t;};vector<C> first;auto all=paths(t);
 for(auto&p:all){T u=native(t,p);first.push_back({dp(u).h.back(),p,u});}
 sort(first.begin(),first.end(),[](auto&a,auto&b){return a.h!=b.h?a.h<b.h:a.p<b.p;});
 cerr<<"source "<<h<<" clone "<<h2<<" paths "<<all.size()<<" minimum "<<first[0].h<<'\n';
 bool found[2]={false,false};I examined=0,hits=0;
 auto frame=[&](const T&p){int i=find(p.begin(),p.end(),port)-p.begin();return pair<int,int>{i?p[i-1]:-1,i+1<n?p[i+1]:-1};};
 for(int rank=start;rank<min(limit,(int)first.size());rank++){
  auto&a=first[rank];cerr<<"rank "<<rank<<" H "<<a.h<<'\n';
  for(auto&q:paths(a.t)){
   if(frame(a.p)!=frame(q))continue;
   T u=native(a.t,q);if(dp(u).h.back()>3*h+4)continue;examined++;
   auto counts=ears(u);
   for(int mask=0;mask<(1<<n);mask++)if(bool(mask>>port&1)==(widened?!incoming:incoming)){
    int lane=counts[mask]==3*h+2?0:counts[mask]==3*h+4?1:-1;
    if(lane<0||found[lane])continue;hits++;T w=birth(u,mask);if(!strong(w))continue;
    int guardian=-1;
    if(widened){
     vector<pair<int,int>> triangles;
     for(int a=0;a<=n;a++)if(w[port]>>a&1)for(int b=0;b<=n;b++)if((w[a]>>b&1)&&(w[b]>>port&1))triangles.push_back({a,b});
     if(triangles.empty())continue;
     for(int candidate:{triangles[0].first,triangles[0].second}){
      bool ok=true;for(auto [a,b]:triangles)if(candidate!=a&&candidate!=b)ok=false;if(ok){guardian=candidate;break;}
     }
     if(guardian<0)continue;
    }
    I hh=dp(clone(w,port)).h.back();if(hh!=3*h2+2+2*lane)continue;
    cout<<"{\"source_H\":"<<h<<",\"port\":"<<port<<",\"offset\":"<<2+2*lane<<",\"clone_source_H\":"<<h2<<",\"clone_target_H\":"<<hh<<",\"rank\":"<<rank<<",\"examined\":"<<examined<<",\"guardian\":"<<guardian<<",\"widened\":"<<(widened?"true":"false")<<",\"first_H\":"<<a.h<<",\"second_H\":"<<dp(u).h.back()<<",\"first_path\":";print(a.p);cout<<",\"second_path\":";print(q);cout<<",\"mask\":"<<mask<<",\"target\":";print(w);cout<<"}\n";cout.flush();found[lane]=true;
    if(found[0]&&found[1])return 0;
   }
  }
 }
 cerr<<"bounded end examined "<<examined<<" scalar hits "<<hits<<" found "<<found[0]<<' '<<found[1]<<'\n';return 0;
}
