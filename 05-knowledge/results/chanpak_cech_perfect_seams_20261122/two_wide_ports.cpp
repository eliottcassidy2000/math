#include <array>
#include <vector>
#include <iostream>
using namespace std;
using Int=long long;
using U=unsigned;
int n,m,buf[12],pc[4]; U g[9],edges[9],ordinary,full; bool ok;
vector<array<int,2>> pairlist;
bool check(int z){
 U out=g[z],in=full^(1<<z)^out;
 vector<U> sh;
 for(U ss=out;ss;ss&=ss-1){
   int s=__builtin_ctz(ss);
   for(U tt=g[s]&in;tt;tt&=tt-1){
     int t=__builtin_ctz(tt);
     U e=((1<<s)|(1<<t))&ordinary;
     if(!e) return false;
     for(U f:sh)if(!(e&f))return false;
     sh.push_back(e);
   }
 }
 return true;
}
void dfs(U todo,int a,int b,int s){
 if(!todo&&a>=1&&b>=1)pc[2*(a-1)+(b-1)]++;
 U next=(todo | (a<2 ? 1u : 0) | (b<2 ? 2u : 0)) & s;
 while(next){int v=__builtin_ctz(next);next&=next-1;dfs(todo&~(1<<v),a+(v==0),b+(v==1),g[v]);}
}
array<int,4> eval(){for(int &x:pc)x=0;dfs(ordinary,0,0,full);return {pc[0]-2*pc[1]-2*pc[2]+4*pc[3],pc[2]-2*pc[3],pc[1]-2*pc[3],pc[3]};}
int main(int argc,char**argv){
 n=atoi(argv[1]);full=(1<<n)-1;ordinary=full^3;
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)pairlist.push_back({i,j});
 // swap marked ports (up to relabelling) to fix 0 -> 1.
 long ct=0,cc=0,rr=0;
 U end=1u<<pairlist.size();
 for(U x=1;x<end;x+=2){
   for(int i=0;i<n;i++)g[i]=0;
   for(unsigned k=0;k<pairlist.size();k++){auto [a,b]=pairlist[k];if(x>>k&1)g[a]|=1<<b;else g[b]|=1<<a;}
   if(!check(0)||!check(1))continue;ct++;
   auto [c,l,r,d]=eval();if(c!=-1||!d||(d&(d-1)))continue;
   vector<array<int,2>> offsets;
   int S=-1-__builtin_ctz(d);
   for(int al=-16;al<=16;al++){int be=S-al;if(be < -20 || be >20)continue;
     if(l*(1LL<<(al+20))+r*(1LL<<(be+20))==-(1LL<<19))offsets.push_back({al,be});
   }
   if(offsets.empty())continue;cc++;
   for(unsigned k=0;k<pairlist.size();k++){
      auto [a,b]=pairlist[k];if(a<2)continue;
      g[a]^=1<<b;g[b]^=1<<a;
      if(check(0)&&check(1)){
        auto [c2,l2,r2,d2]=eval();
        for(auto [al,be]:offsets)if(c2==1&&d2==d&&l2*(1LL<<(al+20))+r2*(1LL<<(be+20))==-(1LL<<19)){
           rr++;cout<<n<<" "<<x<<" flip "<<a<<","<<b<<" offsets "<<al<<","<<be<<" lower "<<c<<","<<l<<","<<r<<","<<d<<" upper "<<c2<<","<<l2<<","<<r2<<","<<d2<<" masks ";
           g[a]^=1<<b;g[b]^=1<<a;for(int i=0;i<n;i++)cout<<g[i]<<",";g[a]^=1<<b;g[b]^=1<<a;cout<<endl;return 0;
        }
      }
      g[a]^=1<<b;g[b]^=1<<a;
   }
 }
 cerr<<"screened "<<ct<<" lower "<<cc<<" hits "<<rr<<endl;
}
