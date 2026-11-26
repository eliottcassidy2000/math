#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>
using namespace std;
struct Opt {int y; array<int,3>s;};
struct RNG {
 uint64_t x=57;
 uint64_t next(){x^=x<<13;x^=x>>7;x^=x<<17;return x;}
 int n(int k){return next()%k;}
 double real(){return (next()>>11)*0x1.0p-53;}
} rng;
int main(int argc,char**argv){
 const int limit=argc>1?atoi(argv[1]):300;
 const int moves=argc>2?atoi(argv[2]):200000;
 int k;cin>>k;
 if(k!=53||limit<1||moves<1)return 2;
 vector<vector<Opt>> a(k);vector<int> base(k,-1),ref(k),branchids;
 for(int h=0;h<57;h++)if(h!=0&&h!=1&&h!=28&&h!=56)branchids.push_back(h);
 // Board order is labelling order from three_stars.SLOTS; assert this outside.
 Opt marker{};
 for(int i=0;i<k;i++){
  int cnt;cin>>ref[i]>>cnt;
  if(branchids[i]==13)ref[i]=16;
  if(branchids[i]==33)ref[i]=19;
  for(int j=0;j<cnt;j++){
   Opt o;cin>>o.y>>o.s[0]>>o.s[1]>>o.s[2];
   if(!i){marker=o;continue;}
   bool legal=true;for(int c=0;c<3;c++){if(o.s[c]==marker.s[c])legal=false;o.s[c]+=56*c;}
   if(!legal)continue;
   if(o.y==ref[i])base[i]=a[i].size();
   a[i].push_back(o);
  }
  if(i&&base[i]<0)return 3;
 }
 int record=4;
 vector<int> best=base;
 for(int loop=0;loop<limit;loop++){
  vector<int> p=(loop%4)?best:base;
  array<int,168> counts{};int c=0;
  auto add=[&](const Opt&o){for(int j:o.s)c+=counts[j]++;};
  auto del=[&](const Opt&o){for(int j:o.s)c-=--counts[j];};
  for(int i=1;i<k;i++)add(a[i][p[i]]);
  if(loop)for(int j=0;j<7+loop%32;j++){
   int i=1+rng.n(k-1);del(a[i][p[i]]);p[i]=rng.n(a[i].size());add(a[i][p[i]]);
  }
  for(int z=0;z<=moves;z++){
   if(c<record){
    record=c;best=p;
    int edits=0;
    for(int i=1;i<k;i++)edits+=p[i]!=base[i]&&branchids[i]!=13&&branchids[i]!=33;
    cout<<"score "<<c<<" edits "<<edits<<" restart "<<loop<<" move "<<z<<'\n';
    for(int i=1;i<k;i++)cout<<branchids[i]<<':'<<a[i][p[i]].y<<' ';
    cout<<endl;
    if(c<=2)return 0;
   }
   int i=1+rng.n(k-1),j=rng.n(a[i].size());
   if(j==p[i])continue;
   int old=c;del(a[i][p[i]]);add(a[i][j]);int delta=c-old;
   // Slowly cool but retain a weak ability to leave plateaus.
   const double T=1.05-0.95*(double)z/moves;
   if(delta<=0||rng.real()<exp(-delta/T))p[i]=j;
   else {del(a[i][j]);add(a[i][p[i]]);}
  }
 }
 cerr<<"finished "<<limit<<" restarts; best "<<record<<" (heuristic only)\n";
}
