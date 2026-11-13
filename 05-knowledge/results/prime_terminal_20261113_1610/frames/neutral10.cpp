// A bounded neutral-rewiring probe around the literal H601 Bezout witness.
#define main census_main
#include "rooted9.cpp"
#undef main
Rows base={2,12,9,17,999,719,911,783,559,15};
vector<pair<int,int>> pairs;
long long tested=0,guarded=0,neutral=0;
vector<tuple<Rows,int,int>> hits;
void visit(Rows &a,int start,int left,int distance){
  if(!left){
    tested++;vector<int> roots;
    for(int r=0;r<10;r++)if(two_run(a,r))roots.push_back(r);
    if(roots.empty())return;guarded++;
    if(hp(a)!=601)return;neutral++;
    for(int r:roots)if(hp(clone(a,r))==1249)hits.emplace_back(a,r,distance);
    return;
  }
  for(int j=start;j<=int(pairs.size())-left;j++){
    auto [u,v]=pairs[j];a[u]^=1u<<v;a[v]^=1u<<u;
    visit(a,j+1,left-1,distance);
    a[u]^=1u<<v;a[v]^=1u<<u;
  }
}
int main(){
  for(int u=0;u<10;u++)for(int v=u+1;v<10;v++)pairs.emplace_back(u,v);
  for(int d=0;d<=4;d++)visit(base,0,d,d);
  cout<<"{\"tested\":"<<tested<<",\"two_run_guarded\":"<<guarded
      <<",\"neutral\":"<<neutral<<",\"hits\":[";
  bool first=true;
  for(auto &[a,r,d]:hits){if(!first)cout<<",";first=false;
    cout<<"{\"root\":"<<r<<",\"distance\":"<<d<<",\"rows\":[";
    for(int i=0;i<10;i++){if(i)cout<<",";cout<<a[i];}cout<<"]}";}
  cout<<"]}\n";
}
